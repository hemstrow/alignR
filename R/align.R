
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
#' @param par numeric, default 1. Number of cores to use for the alignments.
#'
#' @author William Hemstrom
#' @export
align_to_reference <- function(RA_fastqs, RB_fastqs = NULL, reference, par = 1){
  #=============sanity checks===================
  msg <- character()

  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  if(!.check_system_install("samtools")){
    msg <- c(msg, "No samtools install located on system path.\n")
  }
  if(!.check_system_install("bwa")){
    msg <- c(msg, "No bwa install located on system path.\n")
  }

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
  if(!is.null(RB_fastqs)){RB_fastqs <- normalizePath(RB_fastqs)}
  
  if(!is.null(RB)){
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
        cmd <- paste0("bash ", script, " ", chunks[[q]]$RA[i], " ", chunks[[q]]$RB[i], " ", reference)
      }
      system(cmd)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)
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
#' @param RA_fastqs character. Vector of filepaths for the RA (read one) files.
#' @param RB_fastqs character or NULL, default NULL. Vector of filepaths for the
#'   RB (read one) files.
#' @param reference character. Filepath for the reference genome to use for the
#'   alignment.
#' @param M numeric. The \code{-M} parameter from STACKS, the "number of
#'   mismatches allowed between stacks within individuals" (see STACKS
#'   documentation).
#' @param n numeric, default \code{M}. The \code{-n} parameter from STACKS, the
#'   "number of mismatches allowed between stacks between individuals" (see
#'   STACKS documentation).
#' @param par numeric, default 1. Number of cores to use for the alignments.
#' @param check_headers logical, default TRUE. If TRUE, will check that fastq
#'   headers are reasonable (end in /1 and /2 and otherwise match between paired
#'   ends) before running STACKS if doing denovo assembly with paired-end reads.
#'   Somewhat computationally and I/O intensive, so disable if confident that
#'   headers are OK.
#' @param stacks_cleanup logical, default TRUE. If TRUE, files other than .bam
#'   files created by stacks will be removed.
#'   
#' @author William Hemstrom
#' @export
align_denovo <- function(RA_fastqs, RB_fastqs = NULL, M, 
                         n = M, par = 1, check_headers = FALSE,
                         stacks_cleanup = TRUE){
  
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
      
      if(interactive()){
        cat("stacks expects that each RA/RB pair has a unique prefix and end in .1 or .2. alignR can make a 'renamed' directory at the same location as the original files and copy files here with appropriate names.\n")
        resp <- ""
        while(!resp %in% c("y", "n")){
          cat("Would you like to copy and rename files? (y/n)\n")
          resp <- readLines(n = 1)
          resp <- tolower(resp)
        }
        if(resp == "n"){
          stop("stacks expects that each RA/RB pair has a unique prefix.\n")
        }
        else{
          new_dir <- paste0(filepaths[1], "/renamed")
          
          if(!dir.exists(new_dir)){
            dir.create(new_dir)
          }
          
          if(length(unique(map$header) != nrow(map))){
            map$header <- paste0("renamed_", 1:nrow(map))
            map$target_A <- paste0(map$header, ".1.", fileext[1])
            map$target_B <- paste0(map$header, ".2.", fileext[2])
          }
          
          for(i in 1:nrow(map)){
            if(RA_ngz[i] == RA_fastqs[i]){
              system(paste0("cp ", RA_fastqs[i], " ", new_dir, "/", map$target_A[i]))
            }
            else{
              system(paste0("cp ", RA_fastqs[i], " ", new_dir, "/", map$target_A[i], ".gz"))
            }
            
            if(RB_ngz[i] == RB_fastqs[i]){
              system(paste0("cp ", RB_fastqs[i], " ", new_dir, "/", map$target_B[i]))
            }
            else{
              system(paste0("cp ", RB_fastqs[i], " ", new_dir, "/", map$target_B[i], ".gz"))
            }
          }
          
          colnames(map)[5:6] <- c("RA_renamed", "RB_renamed")
          
          rename_key <- map[,c("RA", "RB", "header", "RA_renamed", "RB_renamed")]
          write.table(rename_key, paste0(filepaths[1], "/rename_key.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
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
      else{
        stop("stacks expects that each RA/RB pair has a unique prefix and end in .1 or .2.\n")
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
    cat("Checking and fixing any unrecognized headers in .fastq files.\nProgress:\t")
    
    rstrings <- stringi::stri_rand_strings(length(c(RA_fastqs, RB_fastqs)), length(c(RA_fastqs, RB_fastqs))*10)
    bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
    
    while(any(bad_rstrings)){
      rstrings[bad_rstrings] <- stringi::stri_rand_strings(sum(bad_rstrings), length(c(RA_fastqs, RB_fastqs))*10)
      bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
    }
    
    prog <- 1
    for(i in 1:length(RA_fastqs)){
      cat(i, " ")
      repl_headers(file.path(filepaths[1], map$RA[i]), 
                   file.path(filepaths[1], map$RB[i]), 
                   file.path(filepaths[1], rstrings[prog]),
                   file.path(filepaths[1], rstrings[prog + 1]))
      prog <- prog + 2
    }
    cat("Done.\n")
  }

  # zip if not zipped
  fileext <- tools::file_ext(c(map$RA, map$RB))
  all_files <- file.path(filepaths[1], c(map$RA, map$RB))
  
  if(any(fileext != "gz")){
    if(interactive()){
      cat("stacks expects that each file is zipped (.gz).\n")
      resp <- ""
      while(!resp %in% c("y", "n")){
        cat("Would you like to gzip files? (y/n)\n")
        resp <- readLines(n = 1)
        resp <- tolower(resp)
      }
      if(resp == "n"){
        stop("stacks expects that each file is zipped (.gz).\n")
      }
      
      cat("gzipping", sum(fileext != "gz"), "files... ")
      
      system(paste0("gzip ", paste0(all_files[which(fileext != "gz")], collapse = " ")))
      
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
  
  write.table(popmap, "./alignR_popmap", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  #===========run stacks=====================
  script <- .fetch_a_script("denovo_map_edit.pl", "perl")
  
  cmd <- paste0(script, 
                " --samples ", filepaths[1], 
                " -T ", par, 
                " -M ", M, 
                " -n ", n, 
                " -o ", getwd(), 
                " --popmap alignR_popmap",
                " -e ", ifelse(stacks.ver == "general", "stacks", ""))
  cmd <- ifelse(is_single, cmd, paste0(cmd, " --paired"))
  
  system(cmd)
  
  cat("Cleaning up and renaming if needed.\n")
  if(stacks_cleanup){
    rm.files <- list.files(".", "matches\\.tsv.\\gz")
    rm.files <- c(rm.files, list.files(".", "snps\\.tsv\\.gz"))
    rm.files <- c(rm.files, list.files(".", "alleles\\.tsv\\.gz"))
    rm.files <- c(rm.files, list.files(".", "tags\\.tsv\\.gz"))
    rm.files <- c(rm.files, "tsv2bam.log", "denovo_map.log")
    file.remove(rm.files)
  }
  
  # rename files according to the map
  if(exists("rename_key")){
    bamnames <- paste0(rename_key$header,
                       ".matches.bam")
    new_bamnames <- paste0(tools::file_path_sans_ext(gsub("\\.gz", "", rename_key$RA)), ".matches.bam")
    for(i in 1:nrow(rename_key)){
      mv_cmd <- paste0("mv ", bamnames[i], " ", new_bamnames[i])
      system(mv_cmd)
    }
  }
  else{
    new_bamnames <- paste0(tools::file_path_sans_ext(gsub("\\.gz", "", map$RA)), ".matches.bam")
  }
  
  cat("Indexing reads.\n")
  # index
  for(i in 1:length(new_bamnames)){
    system(paste0("samtools index ", new_bamnames[i]))
  }
  
  cat("Complete!\n")
}
