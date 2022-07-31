.fetch_a_script <- function(script_name, dir) system.file(dir, script_name, package = "alignR")

.check_system_install <- function(code) as.logical(Sys.getenv(paste0(code, "_install")))

.paired_length_check <- function(RA, RB = NULL){
  if(!is.null(RB)){
    return(length(RA) == length(RB))
  }
  else{
    return(TRUE)
  }
}

.check_is_genome <- function(reference){
  if(!tools::file_ext(reference) %in% c("fa", "fna", "fasta")){
    return(paste0("File ", reference, "not .fa, .fna, or .fasta. Is this a genome?\n"))
  }
  else{
    return(TRUE)
  }
}

.angsd_genotyper_key <- function(genotyper){
  genotyper_table <- data.frame(matrix(c(1, "SAMtools",
                                         2, "GATK",
                                         3, "SOAPsnp",
                                         4, "SYK",
                                         5, "phys",
                                         6, "sample"),
                                       ncol = 2, byrow = TRUE))
  if(!genotyper %in% genotyper_table[,2]){
    msg <- c(msg, paste0("Genotyper ", genotyper, " not available. Is this a typo?\n"))
    return(msg)
  }
  else{
    genotyper <- as.numeric(genotyper_table[match(genotyper, genotyper_table[,2]),1])
    return(genotyper)
  }
}

.filter_paralogs <- function(bamfiles, reference,
                             populations = rep("only_pop", length(bamfiles)),
                             outfile = "selected_clean_regions.rf",
                             genotyper = "GATK",
                             SNP_pval = 0.00000001,
                             minQ = 20,
                             minMapQ = 20,
                             buffer = 1000,
                             alpha = pchisq(10, 1, lower.tail = FALSE),
                             par = 1,
                             rf = rf,
                             cleanup = TRUE){
  browser()
  
  #==========setup==============
  pops <- unique(populations)
  if(par > length(pops)){
    if(length(pops) > 1){
      par <- length(pops)
    }
  }
  
  bamfiles <- normalizePath(bamfiles)
  
  #=========prep bamlists==============
  bamlists <- paste0(pops, "_bamlist")
  for(i in 1:length(pops)){
    write.table(bamfiles[which(populations == pops[i])], bamlists[i], quote = F, row.names = FALSE, col.names = FALSE)
  }
  
  #========single pop function============================
  run_one_pop_check <- function(bamfile, genotyper, SNP_pval, minQ, minMapQ, par = 1){
    script <- .fetch_a_script("check_paralogs_single_pop.sh", "shell")
    cmd <- paste0("bash ", script, " ",
                  bamfile, " ",
                  reference, " ",
                  genotyper, " ",
                  SNP_pval, " ",
                  minMapQ, " ",
                  minQ, " ",
                  par)
    
    system(cmd)
  }
  
  #========run check paralogs on each bamlist=============
  if(length(pops) > 1){
    # register parallel architecture
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)
    
    # divide up into ncore chunks
    iters <- length(pops)
    it_par <- (1:iters)%%par
    chunks <- split(bamlists, it_par)
    
    # run
    output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {
      for(i in 1:length(chunks[[q]])){
        run_one_pop_check(chunks[[q]][i], genotyper, SNP_pval, minQ, minMapQ, par = 1)
      }
    }
    
    # release cores and clean up
    parallel::stopCluster(cl)
  }
  else{
    run_one_pop_check(bamlists[1], genotyper, SNP_pval, minQ, minMapQ, par = par)
  }
  
  #========condense=======================================
  # chr lengths
  script <- .fetch_a_script("get_chr_lengths.sh", "shell")
  system(paste0("bash ", script, " ", reference))
  
  # find good regions
  paralog_result_files <- paste0("results_paralogs_", bamlists)
  g.list <- .find_non_paralogous_sections(flist = paralog_result_files, 
                                          unique_chromsome_file = "chr_lengths.txt",
                                          base_pair_buffer = buffer, 
                                          alpha = alpha, existing_rf = rf)
  
  # save
  outfile <- file.path(dirname(bamfiles[1]), outfile)
  write.table(data.frame(p = g.list), file = outfile, col.names = F, row.names = F, quote = F)
  
  # cleanup
  if(cleanup){
    rm.files <- c(paralog_result_files, bamlists,
                  "chr_lengths.txt",
                  paste0("results_snps_", bamlists, ".arg"),
                  paste0("results_snps_", bamlists, ".mafs"),
                  paste0("results_snps_", bamlists, ".pos"),
                  paste0("results_depth_", bamlists))
    file.remove(rm.files)
  }
  
  return(outfile)
}

.find_non_paralogous_sections <- function(flist, unique_chromsome_file, base_pair_buffer, alpha, existing_rf = FALSE){
  browser()

  dat <- vector("list", length(flist))
  for(i in 1:length(dat)){
    dat[[i]] <- read.table(flist[i], header = FALSE)
  }
  dat <- dplyr::bind_rows(dat)

  
  #====================prepare paralog list=====================
  # filter to likely paralogs
  dat[,5] <- pchisq(dat[,5], 1, lower.tail = FALSE)
  
  paralog_list <- dat[which(dat[,5] <= alpha), 1:2]
  rm(dat)
  # unique entries
  paralog_list <- unique(paralog_list)
  colnames(paralog_list) <- c("chr", "position")
  paralog_list <- dplyr::arrange(paralog_list, chr, position) # sort
  
  # now need to highlight the regions to run. These exclude each site within 1kb of a paralog (basically the tag)
  # to do so, loop through each possible chr and find the acceptable regions
  chr_opts <- read.table(unique_chromsome_file, header = F, stringsAsFactors = F)
  chr_lengths <- chr_opts[,2]
  chr_opts <- chr_opts[,1]
  chr_opts <- gsub(">", "", chr_opts)
  
  # if an existing rf is provided, figure out which chrs are actually useable
  if(!isFALSE(existing_rf)){
    rf <- readLines(existing_rf)
    rf_chrs <- gsub(":.+$", "", rf)
    rf_chrs <- gsub(":$", "", rf)
    rf_pos <- gsub("^.+:", "", rf)
    rf_start <- gsub("-.+$", "", rf_pos)
    rf_end <- gsub("^.+-", "", rf_pos)
    rf_start <- gsub("-", "", rf_start)
    rf_end <- gsub("-", "", rf_end)
    chr_opts <- chr_opts[which(chr_opts %in% rf_chrs)] # no reason to even look at any chrs that are entirely excluded based on rf file.
    
    rf <- data.frame(chr = rf_chrs, start = rf_start, end = rf_end)
    
  }
  
  # loop through each possible chromosome, building a series of arguments listing "ok" sites
  good.sections <- character(0)
  for(i in 1:length(chr_opts)){
    
    # if this chr is all good, note and skip
    if(!chr_opts[i] %in% paralog_list$chr){
      good.sections <- c(good.sections, paste0(chr_opts[i], ":1-"))
      next
    }
    
    # initialize
    t.chr.bads <- unlist(paralog_list[which(paralog_list$chr == chr_opts[i]), 2])
    t.chr.bads <- sort(t.chr.bads)
    b.end <- 0
    g.start <- 1
    
    # for each bad entry, need to set a 1kb restriction area
    for(j in 1:length(t.chr.bads)){
      
      # if this is the last one, everything after this bad window is good
      if(j == length(t.chr.bads)){
        # as long as there are actually more bps in the chr, note that everything else is good
        if(chr_lengths[i] >= t.chr.bads[j] + base_pair_buffer + 1){
          else{
            good.sections <- c(good.sections, paste0(chr_opts[i], ":", t.chr.bads[j] + base_pair_buffer + 1, "-"))
          }
        }
        next
      }
      
      # if this bad window start is greater than the previous bad window end, time to print a new good section (unless there is a perfect overlap!)
      if(t.chr.bads[j] - base_pair_buffer > b.end){
        if(g.start < t.chr.bads[j] - base_pair_buffer - 1){
          good.sections <- c(good.sections, paste0(chr_opts[i], ":", g.start, "-", t.chr.bads[j] - base_pair_buffer - 1))
        }
        g.start <- t.chr.bads[j] + base_pair_buffer + 1
        b.end <- t.chr.bads[j] + base_pair_buffer
      }
      
      # otherwise we are still in the same window
      else{
        b.end <- t.chr.bads[j] + base_pair_buffer
        g.start <- t.chr.bads[j] + base_pair_buffer + 1
        next()
      }
    }
  }
  
  browser()
  # harmonize with the existing RF
  
  # remove windows outside of existing rf regions
  
  # truncate windows that are partially inside existing rf regions
  
  # leave other windows intact
  
  return(good.sections)
}
