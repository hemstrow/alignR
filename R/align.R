
#' Filter and align single or paired-end reads to a reference genome.
#'
#' Aligns split, single- or paired-end reads to a reference genome using the bwa
#' mem algorithm, them filters the resulting alignments using a few different
#' approaches via samtools.
#'
#' Paired-end reads are filtered following alignment by first removing PCR
#' duplicates using \code{samtools fixmate} and \code{samtools markdup}, then
#' poorly mapped reads using \code{samtools view} with \code{-q 5}, and then
#' improperly paired reads using \code{samtools view} with  \code{-f 0x2}. The
#' full alignment script can be accessed using \code{system.file("shell",
#' "run_align.sh", "alignR")}.
#'
#' Single-end reads are not filtered following alignment, although filtering can
#' (and usually should) be done while calling genotypes or doing other
#' down-stream analyses. The full alignment script can be accessed using
#' \code{system.file("shell", "run_align_single.sh", "alignR")}.
#'
#' @param RA_fastqs character. Vector of filepaths for the RA (read one) files.
#' @param RB_fastqs character or NULL, default NULL. Vector of filepaths for the
#'   RB (read one) files.
#' @param reference character. Filepath for the reference genome to use for the
#'   alignment.
#' @param mapQ numeric, default 5. Reads with mapping qualities lower than mapQ
#'   will be removed using \code{samtools view -q mapQ}. Only valid for
#'   paired-end reads.
#' @param remove_duplicates logical, default TRUE. If TRUE, PCR duplicates
#'   (clones) will be removed using \code{samtools markdup}. Only valid for
#'   paired-end reads.
#' @param remove_improper_pairs logical, default TRUE. If TRUE, improperly
#'   paired reads reads will be removed using \code{samtools view -f 0x2}. Only
#'   valid for paired-end reads.
#' @param par numeric, default 1. Number of cores to use for the alignments.
#' 
#' @return Generates 'x.sort.flt.bam' and 'x.sort.flt.bam.bai' files (sorted
#'   alignments and index files, respectively), where 'x' is the prefix for the
#'   RA. fastqs. In R, returns a vector of output file paths.
#'   
#' @author William Hemstrom
#' @export
align_reference <- function(RA_fastqs, RB_fastqs = NULL, reference, mapQ = 5,
                            remove_duplicates = TRUE, remove_improper_pairs = TRUE,
                            par = 1){
  #=============sanity checks===================
  msg <- character()

  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  if(!.check_system_install("samtools")){
    msg <- c(msg, "No SAMtools install located on system path.\n")
  }
  if(!.check_system_install("bwa")){
    msg <- c(msg, "No bwa install located on system path.\n")
  }
  
  
  if(!.paired_length_check(RA_fastqs, RB_fastqs)){
    msg <- c(msg, "RA_fastqs and RB_fastqs must be of equal length.\n")
  }
  
  check <- .check_is_genome(reference)
  if(is.character(check)){msg <- c(msg, check)}
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  # check that the reference exists and is indexed.
  if(file.exists(reference)){
    index_files <- paste0(reference, c(".amb", ".ann", ".bwt", ".pac",".sa"))
    if(any(!file.exists(index_files))){
      cat("Genome not indexed. Indexing with 'bwa index'.\n")
      system(paste0("bwa index ", reference))
    }
  }
  else{
    msg <- c(msg, "Reference genome file not found.\n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #==============prepare to run=================
  reference <- normalizePath(reference)
  RA_fastqs <- normalizePath(RA_fastqs)
  RA_gz_files <- which(tools::file_ext(RA_fastqs) == "gz")
  if(length(RA_gz_files) > 0){
    system(paste0("gunzip ", paste0(RA_fastqs[RA_gz_files], collapse = " ")))
    RA_fastqs[RA_gz_files] <- gsub("\\.gz", "", RA_fastqs[RA_gz_files])
  }
  if(!is.null(RB_fastqs)){
    RB_fastqs <- normalizePath(RB_fastqs)
    RB_gz_files <- which(tools::file_ext(RB_fastqs) == "gz")
    if(length(RB_gz_files) > 0){
      system(paste0("gunzip ", paste0(RB_fastqs[RB_gz_files], collapse = " ")))
      RB_fastqs[RB_gz_files] <- gsub("\\.gz", "", RB_fastqs[RB_gz_files])
    }
  }
  
  
  
  if(!is.null(RB_fastqs)){
    is_single <- FALSE
    script <- .fetch_a_script("run_align.sh", "shell")
    
    # prepare file handles
    fastqs <- data.frame(RA = RA_fastqs, RB = RB_fastqs)
    fastqs$RA <- gsub("\\.fastq$", "", fastqs$RA)
    fastqs$RB <- gsub("\\.fastq$", "", fastqs$RB)
  }
  else{
    is_single <- TRUE
    script <- .fetch_a_script("run_align_single.sh", "shell")
    
    # prepare file handles
    fastqs <- data.frame(RA = RA_fastqs)
    fastqs$RA <- gsub("\\.fastq$", "", fastqs$RA)
  }
  

  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # divide up into ncore chunks
  iters <- length(RA_fastqs)
  it_par <- (1:iters)%%par
  chunks <- split(fastqs, it_par)

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:nrow(chunks[[q]])){
      if(is_single){
        cmd <- paste0("bash ", script, " ", chunks[[q]]$RA[i], " ", reference)
      }
      else{
        cmd <- paste0("bash ", script, " ", chunks[[q]]$RA[i], " ", 
                      chunks[[q]]$RB[i], " ", reference, " ",
                      mapQ, " ",
                      ifelse(remove_duplicates, 1, 0), " ",
                      ifelse(remove_improper_pairs, 1, 0))
      }
      system(cmd)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)
  
  return_files <- paste0(fastqs$RA, ".sort.flt.bam")
}


#' Filter and align single or paired-end to a denovo assembly via stacks.
#'
#' Aligns split, single- or paired-end reads to a reference genome using the
#' stacks toolkit. This script is a wrapper for a slightly edited
#' \code{denovo_map.pl} that is adjusted to not use \code{gstacks} or
#' \code{populations}, the former of which can be run with
#' \code{\link{genotype_bams}}. The functionality of the later of which can be
#' done using several other R packages.
#' 
#' Paired-end reads are filtered following alignment by first removing PCR
#' duplicates using \code{samtools fixmate} and \code{samtools markdup}, then
#' poorly mapped reads using \code{samtools view} with \code{-q 5}, and then
#' improperly paired reads using \code{samtools view} with  \code{-f 0x2}. The
#' full alignment script can be accessed using \code{system.file("shell",
#' "run_align.sh", "alignR")}.
#'
#' Single-end reads are not filtered following alignment, although filtering can
#' (and usually should) be done while calling genotypes or doing other
#' down-stream analyses. The full alignment script can be accessed using
#' \code{system.file("shell", "run_align_single.sh", "alignR")}.
#' 
#' Selecting the corretc M and n values for STACKS can be quite tricky.
#' For advice on this, try checking out
#' [this paper by Paris et al 2017](https://doi.org/10.1111/2041-210X.12775),
#' which provides an excellent summary. Make sure to cite this paper if you
#' use it!
#'
#' @param RA_fastqs character. Vector of filepaths for the RA (read one) files.
#' @param RB_fastqs character or NULL, default NULL. Vector of filepaths for the
#'   RB (read one) files.
#' @param M numeric. The \code{-M} parameter from STACKS, the "number of
#'   mismatches allowed between stacks within individuals" (see STACKS
#'   documentation).
#' @param n numeric, default \code{M}. The \code{-n} parameter from STACKS, the
#'   "number of mismatches allowed between stacks between individuals" (see
#'   STACKS documentation).
#' @param par numeric, default 1. Number of cores to use for the alignments.
#' @param check_headers logical, default FALSE. If TRUE, will check that fastq
#'   headers are reasonable (end in /1 and /2 and otherwise match between paired
#'   ends) before running STACKS if doing denovo assembly with paired-end reads.
#'   Somewhat computationally and I/O intensive, so disable if you are confident
#'   that your headers are OK. Headers produced by \code{link{demultiplex}} do
#'   not need checked unless it's \code{stacks_header} argument was set to
#'   FALSE when run.
#' @param stacks_cleanup logical, default TRUE. If TRUE, files other than .bam
#'   files created by stacks will be removed.
#' @param re_align logical, default TRUE. If TRUE, bam files will be re-aligned
#'   to the denovo reference generated by STACKS using \code{bwa mem} via
#'   \code{align_reference}. This is useful if "truncated file" errors are
#'   generated during indexing, which can occasionally occur when stacks aligns
#'   reads.
#' @param mapQ numeric, default 5. Used only if re_align = TRUE, passed to
#'   \code{align_reference}. Reads with mapping qualities lower than mapQ will
#'   be removed using \code{samtools view -q mapQ}. Only valid for paired-end
#'   reads.
#' @param remove_duplicates logical, default TRUE. Used only if re_align = TRUE,
#'   passed to \code{align_reference}. If TRUE, PCR duplicates (clones) will be
#' removed using \code{samtools markdup}. Only valid for paired-end reads.
#' @param remove_improper_pairs logical, default FALSE. Used only if re_align =
#'   TRUE, passed to \code{align_reference}. If TRUE, improperly paired reads
#'   reads will be removed using \code{samtools view -f 0x2}. Only valid for
#'   paired-end reads, and should usually be FALSE for denovo alignments.
#' @param ask_confirmation logical, default TRUE. If TRUE, \code{alignR} will
#'   ask for confirmation before moving or renaming files. If FALSE,
#'   \code{alignR} will move or rename without confirmation.
#'
#' @return If the \code{stacks_cleanup} arugment is TRUE, catalog.fa and index
#'   files for the denovo genome in the directory containing the fastq files
#'   alongside alignments, either 'x.alns.bam', 'x.alns.sort.bam', and
#'   'x.alns.sort.bam.bai' or 'x.sort.flt.bam' and 'x.sort.flt.bam.bai' where
#'   'x' is the prefix for the RA fastqs, depending on if re-alignment with
#'   \code{bwa mem} via \code{align_reference} is requested with the
#'   \code{re_align} argument. In R, returns a vector of output file paths.
#'   
#' @author William Hemstrom
#' @export
align_denovo <- function(RA_fastqs, RB_fastqs = NULL, M, 
                         n = M, par = 1, check_headers = FALSE,
                         stacks_cleanup = TRUE,
                         re_align = TRUE,
                         mapQ = 5,
                         remove_duplicates = TRUE, 
                         remove_improper_pairs = FALSE,
                         ask_confirmation = TRUE){

  #==========sanity checks==========================
  RA_fastqs <- normalizePath(RA_fastqs)
  if(!is.null(RB_fastqs)){RB_fastqs <- normalizePath(RB_fastqs)}
  
  RA_ngz <- gsub("\\.gz$", "", RA_fastqs)
  RB_ngz <- gsub("\\.gz$", "", RB_fastqs)
  
  msg <- character(0)
  
  if(any(duplicated(c(RA_ngz, RB_ngz)))){
    msg <- c(msg, "Some likely duplicated file names provided. Is it possible that you passed both .(fastq).gz and .(fastq) versions of the same file to RA and/or RB?\n")
  }
  
  
  
  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  if(!.check_system_install("samtools")){
    msg <- c(msg, "No samtools install located on system path.\n")
  }
  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }
  stacks.ver <- Sys.getenv("stacks_install")
  if(isFALSE(stacks.ver)){
    msg <- c(msg, "No stacks install located on system path.\n")
  }
  
  if(re_align){
    if(!.check_system_install("bwa")){
      msg <- c(msg, "No bwa install located on system path. This is needed if `re_align` is TRUE.\n")
    }
  }
  
  if(!.paired_length_check(RA_fastqs, RB_fastqs)){
    msg <- c(msg, "RA_fastqs and RB_fastqs must be of equal length.\n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  #===========make popmap and prep===================
  if(!is.null(RB_fastqs)){
    is_single <- FALSE
    
    # prepare file handles
    fastqs <- data.frame(RA = RA_ngz, RB = RB_ngz)
    filepaths <- dirname(c(RA_fastqs, RB_fastqs))
    if(length(unique(filepaths)) != 1){
      stop("stacks expects that all .fastq files are in the same directory. Please move files as needed.\n")
    }
    fileext <- tools::file_ext(c(RA_ngz, RB_ngz))
    if(length(unique(fileext)) != 1){
      stop("stacks expects that all files have the same extension/are of the same type.\n")
    }
    if(!fileext[1] %in% c("fq", "fastq", "fa", "fasta")){
      stop("stack expects that all input files are either .fq, .fastq, .fa, or .fasta.\n")
    }

    map <- data.frame(RA = basename(fastqs$RA),
                      RB = basename(fastqs$RB))
    map$header <- ""
    map$good <- FALSE
    map$target_A <- ""
    map$target_B <- ""
    for(i in 1:nrow(map)){
      key1 <- unlist(strsplit(map$RA[i], ""))
      key2 <- unlist(strsplit(map$RB[i], ""))
      match_key <- key1 == key2
      map$header <- paste0(key1[1:(min(which(!match_key)) - 1)], collapse = "")
      map$target_A <- paste0(map$header[i], ".1.", fileext[i])
      map$target_B <- paste0(map$header[i], ".1.", fileext[i])
    }
    
    
    if(length(unique(map$header)) != nrow(map) | any(map$target_A != RA_ngz) | any(map$target_B != RB_ngz)){
      
      # ask for confirmation if 
      if(interactive() & ask_confirmation){
        cat("stacks expects that each RA/RB pair has a unique prefix and end in .1 or .2. alignR can make a 'renamed' directory at the same location as the original files and copy files here with appropriate names.\n")
        resp <- ""
        while(!resp %in% c("y", "n")){
          if(ask_confirmation){
            cat("Would you like to copy and rename files? (y/n)\n")
            resp <- readLines(n = 1)
            resp <- tolower(resp)
          }
          else{
            resp <- "y"
          }
        }
      }
      else if(!interactive() & ask_confirmation){resp <- "n"}
      else{
        resp <- "y"
      }
      
      
      if(resp == "n"){
        stop("stacks expects that each RA/RB pair has a unique prefix and end in .1 or .2.\n")
      }
      else{
        cat("Copying files.\n")
        new_dir <- paste0(filepaths[1], "/renamed")
        
        if(!dir.exists(new_dir)){
          dir.create(new_dir)
        }
        
        if(length(unique(map$header) != nrow(map))){
          map$header <- paste0("renamed_", 1:nrow(map))
          map$target_A <- paste0(map$header, ".1.", fileext[1])
          map$target_B <- paste0(map$header, ".2.", fileext[2])
        }
        
        
        cl <- parallel::makePSOCKcluster(par)
        doParallel::registerDoParallel(cl)
        output <- foreach::foreach(q = 1:nrow(map), .inorder = TRUE, .errorhandling = "pass",
                                   .packages = "alignR"
        ) %dopar% {
          if(RA_ngz[q] == RA_fastqs[q]){
            system(paste0("cp ", RA_fastqs[q], " ", new_dir, "/", map$target_A[q]))
          }
          else{
            system(paste0("cp ", RA_fastqs[q], " ", new_dir, "/", map$target_A[q], ".gz"))
          }
          
          if(RB_ngz[q] == RB_fastqs[q]){
            system(paste0("cp ", RB_fastqs[q], " ", new_dir, "/", map$target_B[q]))
          }
          else{
            system(paste0("cp ", RB_fastqs[q], " ", new_dir, "/", map$target_B[q], ".gz"))
          }
        }
        
        parallel::stopCluster(cl)
        
        
        colnames(map)[5:6] <- c("RA_renamed", "RB_renamed")
        
        rename_key <- map[,c("RA", "RB", "header", "RA_renamed", "RB_renamed")]
        utils::write.table(rename_key, paste0(filepaths[1], "/rename_key.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
        cat("Renaming successfull. Key for renamed files located at", paste0(filepaths[1], "/rename_key.txt"), "\n\tNote that resulting bam files will be renamed back to their original RA prefixes!\n")
        
        filepaths <- rep(new_dir, length(filepaths))
        map$RA <- map$RA_renamed
        map$RB <- map$RB_renamed
        
        if(any(RA_ngz != RA_fastqs)){
          map$RA[which(RA_ngz != RA_fastqs)] <- paste0(map$RA[which(RA_ngz != RA_fastqs)], ".gz")
        }
        if(any(RB_ngz != RB_fastqs)){
          map$RB[which(RB_ngz != RB_fastqs)] <- paste0(map$RB[which(RB_ngz != RB_fastqs)], ".gz")
        }
      }
    }
    
    popmap <- cbind(header = map$header, pop = "filler")
  }
  else{
    popmap <- cbind(header = RA_ngz, pop = "filler")
    map <- data.frame(RA = RA_fastqs)
  }
  
  
  # check for bad headers
  repl_headers <- function(file1, file2, rstring, rstring2){
    script <- .fetch_a_script("header_fixer_paired.pl", "perl")
    
    
    if(tools::file_ext(file1) == "gz"){
      rezip1 <- TRUE
      system(paste0("gunzip ", file1))
      file1 <- tools::file_path_sans_ext(file1)
    }
    else{
      rezip1 <- FALSE
    }
    
    if(tools::file_ext(file2) == "gz"){
      rezip2 <- TRUE
      system(paste0("gunzip ", file2))
      file2 <- tools::file_path_sans_ext(file2)
    }
    else{
      rezip2 <- FALSE
    }
    
    system(paste0("perl ", script, " ", file1, " ", file2, " ", rstring, " ", rstring2))
    system(paste0("mv ", rstring, " ", file1))
    system(paste0("mv ", rstring2, " ", file2))
    
    if(rezip1 | rezip2){
      system(paste0("gzip ", ifelse(rezip1, file1,""), " ", ifelse(rezip2, file2, "")))
    }
  }
  
  if(check_headers & !is.null(RB_fastqs)){
    cat("Checking and fixing any unrecognized headers in .fastq files.\n")
    
    rstrings <- .rand_strings(length(c(RA_fastqs, RB_fastqs)), 10)
    bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
    
    while(any(bad_rstrings)){
      rstrings[bad_rstrings] <- .rand_strings(sum(bad_rstrings), length(c(RA_fastqs, RB_fastqs)))
      bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
    }
    
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)
    output <- foreach::foreach(q = 1:nrow(map), .inorder = TRUE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {
      repl_headers(file.path(filepaths[1], map$RA[q]), 
                   file.path(filepaths[1], map$RB[q]), 
                   file.path(filepaths[1], rstrings[q]),
                   file.path(filepaths[1], rstrings[(length(rstrings)/2) + q]))
      TRUE
    }
    parallel::stopCluster(cl)
    
    cat("Done.\n")
  }

  # zip if not zipped
  fileext <- tools::file_ext(c(map$RA, map$RB))
  all_files <- file.path(filepaths[1], c(map$RA, map$RB))
  
  if(any(fileext != "gz")){
    if(interactive() & ask_confirmation){
      cat("stacks expects that each file is zipped (.gz).\n")
      resp <- ""
      while(!resp %in% c("y", "n")){
        if(ask_confirmation){
          cat("Would you like to gzip files? (y/n)\n")
          resp <- readLines(n = 1)
          resp <- tolower(resp)
        }
        else{
          resp <- "y"
        }
      }
    }
    else if(!interactive() & ask_confirmation){
      resp <- "n"
    }
    else{
      resp <- "y"
    }
    
    if(resp == "y"){
      
      cat("gzipping", sum(fileext != "gz"), "files... ")
      
      cl <- parallel::makePSOCKcluster(par)
      doParallel::registerDoParallel(cl)
      output <- foreach::foreach(q = 1:sum(fileext != "gz"), .inorder = TRUE, .errorhandling = "pass",
                                 .packages = "alignR"
      ) %dopar% {
        system(paste0("gzip ", all_files[which(fileext != "gz")][q]))
      }
      parallel::stopCluster(cl)
      cat("Done.\n")
      
      map$RA[which(tools::file_ext(map$RA) != "gz")] <-
        paste0(map$RA[which(tools::file_ext(map$RA) != "gz")], ".gz")
      
      if(!is.null(RB_fastqs)){
        map$RB[which(tools::file_ext(map$RB) != "gz")] <-
          paste0(map$RB[which(tools::file_ext(map$RB) != "gz")], ".gz")
      }
    }
    else{
      stop("stacks expects that each file is zipped (.gz).\n")
    }
  }
  
  utils::write.table(popmap, "./alignR_popmap", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  #===========run stacks=====================
  script <- .fetch_a_script("denovo_map_edit.pl", "perl")
  
  cmd <- paste0(script, 
                " --samples ", filepaths[1], 
                " -T ", par, 
                " -M ", M, 
                " -n ", n, 
                " -o ", filepaths[1], 
                " --popmap alignR_popmap",
                " -e ", ifelse(stacks.ver == "general", "stacks", ""))
  cmd <- ifelse(is_single, cmd, paste0(cmd, " --paired"))
  
  system(cmd)
  gcmd <- paste0(ifelse(stacks.ver == "general", "stacks ", ""),
                 "gstacks -P ", 
                 filepaths[1],
                 " -M alignR_popmap -t ",
                 par,
                 ifelse(re_align, "", " --write-alignments"))
  system(gcmd)
  
  #===========clean-up and post-processing===========
  
  cat("Cleaning up and renaming if needed.\n")
  if(stacks_cleanup){
    rm.files <- list.files(filepaths[1], "matches\\.tsv.\\gz$", full.names = TRUE)
    rm.files <- c(rm.files, list.files(filepaths[1], "matches\\.bam$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "snps\\.tsv\\.gz$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "alleles\\.tsv\\.gz$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "tags\\.tsv\\.gz$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "snps\\.tsv$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "alleles\\.tsv$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "tags\\.tsv$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "matches\\.tsv$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "tsv2bam\\.log$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "denovo_map\\.log$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "gstacks\\.log$", full.names = TRUE))
    rm.files <- c(rm.files, list.files(filepaths[1], "gstacks\\.log\\.distribs$", full.names = TRUE))
    rm.files <- c(rm.files, "alignR_popmap")
    file.remove(rm.files)
  }
  

  if(re_align){
    # index
    cat("Re-aligning reads using `align_reference`.\n")
    system(paste0("gunzip ", file.path(filepaths[1], "catalog.fa.gz")))
    system(paste0("bwa index ", file.path(filepaths[1], "catalog.fa")))
    
    return_files <- align_reference(RA_fastqs = RA_fastqs, 
                                    RB_fastqs = RB_fastqs,
                                    reference = file.path(filepaths[1], "catalog.fa"), 
                                    par = par, 
                                    mapQ = mapQ, 
                                    remove_duplicates = remove_duplicates, 
                                    remove_improper_pairs = remove_improper_pairs)
  }
  
  else{
    # rename files according to the map
    if(exists("rename_key")){
      bamnames <- file.path(filepaths[1], paste0(rename_key$header,
                         ".alns.bam"))
      new_bamnames <- file.path(filepaths[1],
                                paste0(tools::file_path_sans_ext(gsub("\\.gz", "", rename_key$RA)), ".alns.bam"))
      for(i in 1:nrow(rename_key)){
        mv_cmd <- paste0("mv ", bamnames[i], " ", new_bamnames[i])
        system(mv_cmd)
      }
    }
    else{
      new_bamnames <- file.path(filepaths[1],
                                paste0(tools::file_path_sans_ext(gsub("\\.gz", "", map$RA)), ".alns.bam"))
    }
    
    cat("Indexing reads with SAMtools.\n")
    sort_names <- gsub("\\.bam$", ".sort.bam", new_bamnames)
    for(i in 1:length(new_bamnames)){
      system(paste0("samtools sort -o ", sort_names[i], " ", new_bamnames[i]))
      system(paste0("samtools index ", sort_names[i]))
    }
    
    return_files <- sort_names
  }
  


  cat("Complete!\n")
  return(return_files)
}
