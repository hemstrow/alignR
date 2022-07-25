
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
#' @param par numeric, default 1. Number of cores to use for the alignments.
#'
#' @author William Hemstrom
#' @export
align_denovo <- function(RA_fastqs, RB_fastqs = NULL, M, n, par = 1){
  
  #===========make popmap===================
  RA_fastqs <- normalizePath(RA_fastqs)
  if(!is.null(RB_fastqs)){RB_fastqs <- normalizePath(RB_fastqs)}
  
  RA_ngz <- gsub("\\.gz$", "", RA_fastqs)
  RB_ngz <- gsub("\\.gz$", "", RB_fastqs)
  
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
          
          write.table(map[,c("RA", "RB", "header", "RA_renamed", "RB_renamed")], paste0(filepaths[1], "/rename_key.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
          cat("Renaming successfull. Key for renamed files located at", paste0(filepaths[1], "/rename_key.txt"), "\n")
          
          filepaths <- rep(new_dir, length(filepaths))
          map$RA <- map$RA_renamed
          map$RB <- map$RB_renamed
        }
      }
      else{
        stop("stacks expects that each RA/RB pair has a unique prefix and end in .1 or .2.\n")
      }
      
    }
    
    popmap <- cbind(header = map$header, pop = "filler")
  }
  else{

  }
  
  
  # check for bad headers
  repl_headers <- function(file, rstring){
    script <- .fetch_a_script("header_fixer.pl", "perl")
    
    
    if(tools::file_ext(file) == "gz"){
      rezip <- TRUE
      system("gunzip ", file)
    }
    else{
      rezip <- FALSE
    }
    
    system(paste0("perl ", script, " ", file, " ", rstring))
    system(paste0("mv ", rstring, " ", file))
    
    if(rezip){
      system("gzip ", file)
    }
  }
  
  cat("Checking and fixing any unrecognized headers in .fastq files.\n")

  rstrings <- stringi::stri_rand_strings(length(c(RA_fastqs, RB_fastqs)), length(c(RA_fastqs, RB_fastqs))*10)
  bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
  
  while(any(bad_rstrings)){
    rstrings[bad_rstrings] <- stringi::stri_rand_strings(sum(bad_rstrings), length(c(RA_fastqs, RB_fastqs))*10)
    bad_rstrings <- file.exists(file.path(filepaths[1], rstrings))
  }
  
  browser()
  prog <- 1
  for(i in 1:length(RA_fastqs)){
    repl_headers(file.path(filepaths[1], map$RA[i]), file.path(filepaths[1], rstrings[prog]))
    prog <- prog + 1
    
    if(!is.null(RB_fastqs)){
      repl_headers(file.path(filepaths[1], map$RB[i]), file.path(filepaths[1], rstrings[prog]))
      prog <- prog + 1
    }
  }
  

  # zip if not zipped
  fileext <- tools::file_ext(c(RA_fastqs, RB_fastqs))
  
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
      
      for(i in 1:length(fileext)){
        if(fileext[i] != "gz"){
          system(paste0("gzip -c ", c(RA_fastqs, RB_fastqs)[i], " > ", paste0(c(RA_fastqs, RB_fastqs)[i], ".gz")))
        }
      }
      
      RA_fastqs[which(tools::file_ext(RA_fastqs) != "gz")] <-
        paste0(RA_fastqs[which(tools::file_ext(RA_fastqs) != "gz")], ".gz")
      
      if(!is.null(RB_fastqs)){
        RB_fastqs[which(tools::file_ext(RB_fastqs) != "gz")] <-
          paste0(RB_fastqs[which(tools::file_ext(RB_fastqs) != "gz")], ".gz")
      }
    }
    else{
      stop("stacks expects that each file is zipped (.gz).\n")
    }
  }
  
  write.table(popmap, "./alignR_popmap", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  #===========run stacks=====================
  browser()
  script <- .fetch_a_script("denovo_map_edit.pl", "perl")
  
  cmd <- paste0(script, 
                " --samples ", filepaths[1], 
                " -T ", par, 
                " -M ", M, 
                " -n ", n, 
                " -o ", getwd(), 
                " --popmap alignR_popmap")
  cmd <- ifelse(is_single, cmd, paste0(cmd, " --paired"))
  
  system(cmd)
}
