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
                             alpha = stats::pchisq(10, 1, lower.tail = FALSE),
                             par = 1,
                             rf = rf,
                             cleanup = TRUE){
  
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
    utils::write.table(bamfiles[which(populations == pops[i])], bamlists[i], quote = F, row.names = FALSE, col.names = FALSE)
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
  utils::write.table(data.frame(p = g.list), file = outfile, col.names = F, row.names = F, quote = F)
  
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
  
  chr <- position <- NULL

  dat <- vector("list", length(flist))
  for(i in 1:length(dat)){
    dat[[i]] <- utils::read.table(flist[i], header = FALSE)
  }
  dat <- dplyr::bind_rows(dat)

  
  #====================prepare paralog list=====================
  # filter to likely paralogs
  dat[,5] <- stats::pchisq(dat[,5], 1, lower.tail = FALSE)
  
  paralog_list <- dat[which(dat[,5] <= alpha), 1:2]
  rm(dat)
  # unique entries
  paralog_list <- unique(paralog_list)
  colnames(paralog_list) <- c("chr", "position")
  paralog_list <- dplyr::arrange(paralog_list, chr, position) # sort
  
  # now need to highlight the regions to run. These exclude each site within 1kb of a paralog (basically the tag)
  # to do so, loop through each possible chr and find the acceptable regions
  chr_opts <- utils::read.table(unique_chromsome_file, header = F, stringsAsFactors = F)
  chr_lengths <- chr_opts[,2]
  chr_opts <- chr_opts[,1]
  chr_opts <- gsub(">", "", chr_opts)
  
  # if an existing rf is provided, figure out which chrs are actually useable
  if(!isFALSE(existing_rf)){
    rf <- readLines(existing_rf)
    rf_chrs <- gsub(":.+$", "", rf)
    rf_chrs <- gsub(":$", "", rf_chrs)
    rf_pos <- gsub("^.+:", "", rf)
    rf_start <- gsub("-.+$", "", rf_pos)
    rf_end <- gsub("^.+-", "", rf_pos)
    rf_start <- gsub("-", "", rf_start)
    rf_end <- gsub("-", "", rf_end)
    rf <- data.frame(chr = rf_chrs, start = rf_start, end = rf_end)
    
    need_starts <- which(rf$start == "")
    rf$start[need_starts] <- chr_lengths[match(rf$chr[need_starts], chr_opts)]
    need_ends <- which(rf$end == "")
    rf$end[need_ends] <- chr_lengths[match(rf$chr[need_ends], chr_opts)]
    rf$start <- as.numeric(rf$start)
    rf$end <- as.numeric(rf$end)
    
    
    keep_chrs <- which(chr_opts %in% rf_chrs)
    chr_opts <- chr_opts[keep_chrs] # no reason to even look at any chrs that are entirely excluded based on rf file.
    chr_lengths <- chr_lengths[keep_chrs]
    
    
    
  }
  
  #======================loop through each possible chromosome, building a series of arguments listing "ok" sites=================
  good.sections <- list(chr = character(), start = numeric(), end = numeric())
  for(i in 1:length(chr_opts)){
    
    # if this chr is all good, note and skip
    if(!chr_opts[i] %in% paralog_list$chr){
      good.sections$chr <- c(good.sections$chr, chr_opts[i])
      good.sections$start <- c(good.sections$start, 1)
      good.sections$end <- c(good.sections$end, chr_lengths[i])
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
          
          good.sections$chr <- c(good.sections$chr, chr_opts[i])
          good.sections$start <- c(good.sections$start, t.chr.bads[j] + base_pair_buffer + 1)
          good.sections$end <- c(good.sections$end, chr_lengths[i])
        }
        next
      }
      
      # if this bad window start is greater than the previous bad window end, time to print a new good section (unless there is a perfect overlap!)
      if(t.chr.bads[j] - base_pair_buffer > b.end){
        if(g.start < t.chr.bads[j] - base_pair_buffer - 1){
          good.sections$chr <- c(good.sections$chr, chr_opts[i])
          good.sections$start <- c(good.sections$start, g.start)
          good.sections$end <- c(good.sections$end, t.chr.bads[j] - base_pair_buffer - 1)
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
  
  #===================harmonize with the existing RF (if provided) via looping through smaller rf================
  # skip everything else if no existing rf to mesh with
  if(isFALSE(existing_rf)){
    return(paste0(good.sections$chr, ":", good.sections$start, "-", good.sections$end))
  }
  
  cat("Merging accepted regions (only regions ID'd as being both non-paralogs and in the provided .rf file accepted).\n This may take a few minutes.\nProgress:\n")
  good.sections <- as.data.frame(good.sections)
  
  # figure out which rf is smaller (quicker to loop through, since the n*p part is vectorized)
  if(nrow(good.sections) < nrow(rf)){
    rf1 <- good.sections
    rf2 <- rf
  }
  else{
    rf1 <- rf
    rf2 <- good.sections
  }
  rm(rf, good.sections)
  
  # loop through each chr, then each region in rf1 to do the comparisons and save the results.
  prog <- 0
  next_completion <- 5
  write.sections <- character()
  for(i in 1:length(chr_opts)){
    rf2.consider <- which(rf2$chr == chr_opts[i])
    rf1.consider <- which(rf1$chr == chr_opts[i])
    for(j in 1:length(rf1.consider)){
      prog <- prog + 1
      current_completion <- (prog/nrow(rf1)) * 100
      if(current_completion > next_completion){
        cat(paste0("\t", floor(current_completion), "%\n"))
        next_completion <- next_completion + 5
      }
      
      # these logicals are all <=, since = means we include that bp. Had to plot this shit out.
      l1 <- rf1[rf1.consider[j],]$start <= rf2[rf2.consider,]$start # rf1 start before rf2 start?
      l2 <- rf1[rf1.consider[j],]$start <= rf2[rf2.consider,]$end # rf1 start before rf2 end?
      l3 <- rf1[rf1.consider[j],]$end <= rf2[rf2.consider,]$start # rf1 end before rf2 start?
      l4 <- rf1[rf1.consider[j],]$end <= rf2[rf2.consider,]$end # rf1 end before rf2 end?
      
      # if the r1 section is entirely within an rf2 region, keep the whole r1 and move on
      r1_within_r2 <- which(!l1 & l2 & !l3 & l4)
      if(length(r1_within_r2) != 0){
        write.sections <- c(write.sections, paste0(chr_opts[i], ":", 
                                                   rf1[rf1.consider[j],]$start,
                                                   rf1[rf1.consider[j],]$end))
        next # can move on, since if this is ever true none of the others can be, and all of the other rf2 regions are rejected.
      }
      

      # keep whichever r2 sections are completely inside r1
      r2_within_r1 <- which(l1 & l2 & !l3 & !l4)
      if(length(r2_within_r1) > 0){
        write.sections <- c(write.sections, paste0(chr_opts[i], ":",
                                                   rf2[rf2.consider[r2_within_r1],]$start, "-",
                                                   rf2[rf2.consider[r2_within_r1],]$end))
      }
      
      # truncate where r1 overhangs to the left of r2
      r1_left_truncate <- which(l1 & l2 & !l3 & l4)
      if(length(r1_left_truncate) != 0){
        write.sections <- c(write.sections, paste0(chr_opts[i], ":",
                                                   rf2[rf2.consider[r1_left_truncate],]$start, "-",
                                                   rf1[rf1.consider[j],]$end))
      }
      
      # truncate where r1 overhangs to the right of r2
      r1_right_truncate <- which(!l1 & l2 & !l3 & !l4)
      if(length(r1_right_truncate)){
        write.sections <- c(write.sections, paste0(chr_opts[i], ":",
                                                   rf1[rf1.consider[j],]$start, "-",
                                                   rf2[rf2.consider[r1_right_truncate],]$end))
      }
      
      # TTTT and FFFF have no overlap, don't need to write anything
    }
    
    
  }
  
  return(write.sections)
}


.rand_strings <- function(n, length){
  chrs <- sample(c(LETTERS, letters, 1:9), n*length, replace = TRUE)
  chrs <- split(chrs, rep(1:n, length.out = n*length))
  chrs <- unlist(lapply(chrs, function(x) paste0(x, collapse = "")))
  chrs <- as.character(chrs)
  return(chrs)
}



.dependency_function_match <- function(dependancies){
  dep_tab <- list(bash = c("all"),
                  perl = c("align_denovo", "plate_split", "demultiplex", "genotype_bams with doVcf = TRUE and doGeno = 'NN' or 'numeric'."),
                  stacks = c("align_denovo"),
                  bcftools = c("genotype_bams with doVcf and doGeno other than 'NN' or 'numeric'"),
                  ngsParalog = c("genotype_bams with filter_paralogs = TRUE"),
                  angsd = c("genotype_bams"),
                  samtools = c("align_reference", "align_denovo", "genotype_bams"),
                  bwa = c("align_reference", "align_denovo with re_align = TRUE"))
  
  return(dep_tab[names(dep_tab) %in% dependancies])
}


.rename_files <- function(from, to){
  for(i in 1:length(from)){
    cmd <- paste0("mv ", from[i], " ", to[i])
    system(cmd)
  }
}


# takes a two tab file with chromsome names and lengths and a chunck size
# saves .rf files batched into chunk_size bp lengths for genotyping
# returns a vector of the rf file names
.batch_genotyping <- function(chr_info, chunk_size){
  old_scipen <- options("scipen")
  options(scipen = 999)
  chr_info <- read.table(chr_info, sep = "\t")
  chr_info$cumsum <- cumsum(as.numeric(chr_info$V2))
  chr_info$css <- c(0, chr_info$cumsum[-nrow(chr_info)]) + 1
  
  starts <- seq(1, chr_info$cumsum[nrow(chr_info)], by = chunk_size)
  ends <- c(starts[-1] - 1, chr_info$cumsum[nrow(chr_info)])
  
  rfs <- vector("list", length(starts))
  
  for(i in 1:length(starts)){
    tchrs <- chr_info[chr_info$cumsum >= starts[i] & chr_info$css <= ends[i],]
    start_chr_start <- starts[i] - tchrs$css[1] + 1
    end_chr_end <- ends[i] - tchrs$css[nrow(tchrs)] + 1
    if(nrow(tchrs) == 1){
      rfs[[i]] <- data.frame(chr = tchrs$V1, start = start_chr_start, end = end_chr_end)
    }
    else if(nrow(tchrs) == 2){
      rfs[[i]] <- data.frame(chr = tchrs$V1,
                             start = c(start_chr_start, 1),
                             end =  c(tchrs$V2[-nrow(tchrs)], end_chr_end))
    }
    else{
      
      rfs[[i]] <- data.frame(chr = tchrs$V1,
                             start = c(start_chr_start, rep(1, nrow(tchrs) - 1)),
                             end = c(tchrs$V2[-nrow(tchrs)], end_chr_end))
    }
    
    rfs[[i]]$cmd <- paste0(rfs[[i]]$chr, ":", rfs[[i]]$start, "-", rfs[[i]]$end)
    
    writeLines(rfs[[i]]$cmd, paste0("region", i, ".rf"), sep = "\n")
  }
  
  options(scipen = old_scipen$scipen)
  return(paste0("region", 1:length(starts), ".rf"))
}