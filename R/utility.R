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
  
  return(return_files)
}
