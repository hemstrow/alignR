#' Equalize the number of reads across multiple fastq or fastq files.
#'
#' Subset out the first n reads of every fastq/fastq file, where n is either the
#' smallest number of reads across all samples or user supplied.
#'
#' @param RA_fastqs character. Vector of filepaths for the RA (read one) files.
#' @param RB_fastqs character or NULL, default NULL. Vector of filepaths for the
#'   RB (read one) files. If NULL, assumes reads are single-end.
#' @param nreads numeric or "lowest", default "lowest". Number of reads to keep.
#'   If "lowest", will keep n reads, where n is the number of reads in the
#'   smallest .fastq/.fasta file.
#' @param outfile_prefix character. Prefix to be appended to each resulting
#'   .fastq/.fastq file.
#'
#' @export
#' @author William Hemstrom
#' @return Generates files subset to the specified size. In R, returns a list of
#'   RA (and RB if paired-end) file paths to subset files.
equalize_read_counts <- function(RA_fastqs, RB_fastqs = NULL, nreads = "lowest",
                                 outfile_prefix = "subset_"){
  old.scipen <- options("scipen")
  options(scipen = 999)
  #==========sanity checks============
  RA_fastqs <- normalizePath(RA_fastqs)
  if(!is.null(RB_fastqs)){
    RB_fastqs <- normalizePath(RB_fastqs)
  }

  msg <- character(0)

  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }

  exts <- tools::file_ext(c(RA_fastqs, RB_fastqs))
  if(length(unique(exts)) != 1){
    msg <- c(msg, "RA/RB files must all have identical extensions.\n")
  }
  if(any(!exts %in% c("fa", "fq", "fasta", "fastq"))){
    msg <- c(msg, "Only fasta, fa, fastq, or fq files accepted.\n")
  }


  if(!.paired_length_check(RB_fastqs, RA_fastqs)){
    msg <- c(msg, "RA_fastqs and RB_fastqs must be of equal length.\n")
  }


  if(length(msg) != 0){
    stop(msg)
  }

  #=========subset====================
  if(nreads == "lowest"){
    wc <- numeric(length(RA_fastqs))
    for(i in 1:length(wc)){
      wc[i] <- as.numeric(system(paste0("wc -l < ", RA_fastqs[i]), intern = TRUE))
    }

    nreads <- min(wc)
  }
  else{
    if(tools::file_ext(RA_fastqs[1]) %in% c("fastq", "fq")){
      nreads <- nreads * 4 # convert to number of lines to keep
    }
    else if(tools::file_ext(RA_fastqs[1]) %in% c("fasta", "fa")){
      nreads <- nreads * 2 # convert to number of lines to keep
    }
  }


  RAfp <- file.path(dirname(RA_fastqs), paste0(outfile_prefix, basename(RA_fastqs)))
  if(!is.null(RB_fastqs)){
    RBfp <- file.path(dirname(RB_fastqs), paste0(outfile_prefix, basename(RB_fastqs)))
  }


  for(i in 1:length(RA_fastqs)){
    system(paste0("head -n ", nreads, " ", RA_fastqs[i], " > ", RAfp[i]))
    if(!is.null(RB_fastqs)){
      system(paste0("head -n ", nreads, " ", RB_fastqs[i], " > ", RBfp[i]))
    }
  }

  return_files <- list(RA = RAfp, RB = RBfp)
  if(is.null(return_files$RB)){
    return_files$RB <- NULL
  }

  options(scipen = old.scipen$scipen)
  return(return_files)
}

#' Merge fastq files
#'
#' Merges fastq files as a set from each element of a list. Resulting merged
#' files will be in the directory of the first fastq file of the first set to
#' merge.
#'
#' @param file_list list where each element contains a set of fastq files to
#'   merge.
#' @param names character vector, default NULL. A vector of names for each
#'   merged output file. If not provided, files will be named "merged_x", where
#'   x is the name of the first file in each set to merge.
#' @param par numeric, default 1. Number of cores to use for the merges.
#'
#' @export
#' @author William Hemstrom
#' @return Generates merged fastq files. File paths provided as a character
#'   vector returned from this function.
merge_fastqs <- function(file_list, names = NULL, par = 1){
  #==========sanity checks=============
  msg <- character()

  if(!is.null(names)){
    if(is.character(names)){
      if(length(names) != length(file_list)){
        msg <- c(msg, "The number of provided names must equal the number of desired output files.\n")
      }
    }
    else{
      msg <- c(msg, "The provided names must be a character vector.\n")
    }
  }

  if(any(!tools::file_ext(unlist(file_list)) %in% c("fasta", "fastq", "fa", "fq"))){
    msg <- c(msg, "Only fastq, fastq, fa, or fq files accepted.\n")
  }

  bad.files <- !file.exists(unlist(file_list))
  if(any(bad.files)){
    msg <- c(msg, paste0("Some files not located:\n\t", paste0(unlist(file_list)[bad.files], collapse = "\n\t"), "\n"))
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #==========sanity checks=============

  # prep
  par <- min(par, length(file_list))

  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  dir <- normalizePath(dirname(file_list[[1]][1]))

  # run
  output <- foreach::foreach(q = 1:length(file_list), .inorder = TRUE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {

    if(is.null(names)){
      outfile <- file.path(dir, paste0("merged_", basename(file_list[[q]][1])))
    }
    else{
      outfile <- file.path(dir, paste0(names[q], ".fastq"))
    }

    system(paste0("cat ", paste0(file_list[[q]], collapse = " "), " > ", outfile))

    outfile
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(unlist(output))
}



#' Merge bam files
#'
#' Merges bam files as a set from each element of a list using samtools.
#' Resulting merged files will be in the directory of the first bam file of the
#' first set to merge. Bamfiles must be sorted, see details.
#'
#' Bamfiles must be sorted to merge them using samtools. Outputs from both
#' bamfile-generating alignR functions (\code{align_denovo} and
#' \code{align_reference}) are sorted. Otherwise, bamfiles can be sorted using
#' the \code{sort_bamfiles} utility function.
#'
#' @param file_list list where each element contains a set of bam files to
#'   merge. Each bam file must be sorted, see details.
#' @param names character vector, default NULL. A vector of names for each
#'   merged output file. If not provided, files will be named "merged_x", where
#'   x is the name of the first file in each set to merge.
#' @param par numeric, default 1. Number of cores to use for the merges.
#'
#' @export
#' @author William Hemstrom
#' @return Generates merged bam files. File paths provided as a character
#'   vector returned from this function.
merge_bams <- function(file_list, names = NULL, par = 1){
  #==========sanity checks=============
  msg <- character()

  if(!is.null(names)){
    if(is.character(names)){
      if(length(names) != length(file_list)){
        msg <- c(msg, "The number of provided names must equal the number of desired output files.\n")
      }
    }
    else{
      msg <- c(msg, "The provided names must be a character vector.\n")
    }
  }

  if(any(!tools::file_ext(unlist(file_list)) %in% c("bam"))){
    msg <- c(msg, "Only bam files accepted.\n")
  }

  bad.files <- !file.exists(unlist(file_list))
  if(any(bad.files)){
    msg <- c(msg, paste0("Some files not located:\n\t", paste0(unlist(file_list)[bad.files], collapse = "\n\t"), "\n"))
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #==========sanity checks=============

  # prep
  par <- min(par, length(file_list))

  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  dir <- normalizePath(dirname(file_list[[1]][1]))

  # run
  output <- foreach::foreach(q = 1:length(file_list), .inorder = TRUE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {

    if(is.null(names)){
      outfile <- file.path(dir, paste0("merged_", basename(file_list[[q]][1])))
    }
    else{
      outfile <- file.path(dir, paste0(names[q], ".bam"))
    }

    system(paste0("samtools merge -o ", outfile, " ", paste0(file_list[[q]], collapse = " ")))

    outfile
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(unlist(output))
}

#' Sort bam files
#'
#' Sorts bam files using samtools.
#'
#' @param bams A list of bamfiles to be sorted.
#' @param par numeric, default 1. Number of cores to use for sorting.
#'
#' @export
#' @author William Hemstrom
#' @return Generates sorted bam files. File paths provided as a character
#'   vector returned from this function.
sort_bams <- function(bams, par = 1){
  #==========sanity checks=============
  msg <- character()

  if(any(!tools::file_ext(unlist(bams)) %in% c("bam"))){
    msg <- c(msg, "Only bam files accepted.\n")
  }

  bad.files <- !file.exists(unlist(bams))
  if(any(bad.files)){
    msg <- c(msg, paste0("Some files not located:\n\t", paste0(unlist(bams)[bad.files], collapse = "\n\t"), "\n"))
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #==========sanity checks=============

  # prep
  par <- min(par, length(bams))

  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # run
  output <- foreach::foreach(q = 1:length(bams), .inorder = TRUE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    outfile <- gsub("\\.bam$", "\\.sort\\.bam", bams[q])
    system(paste0("samtools sort -o", outfile, " ", bams[q]))
    outfile
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(unlist(output))
}
