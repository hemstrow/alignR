#' Split paired-end .fastq files with reads from different plates.
#'
#' Splits paired-end sequencing .fastq files containing reads from multiple
#' different plates with unique plate indices into different fastq files for
#' each plate index.
#'
#' The R1 and R3 .fastq files should contain the actual forward and reverse
#' reads. The optional R2 .fastq file should contain only the plate indices. As
#' is standard, the reads are assumed to be in the same order for each of the
#' .fastq files. If no R2 file is provided, the index is assumed to be at the
#' end of the headers for each read in the R1 and R3 files.
#'
#' @param R1 character. File name for R1 (read one) file.
#' @param R2 character or NULL, default NULL. File name for R2 (index read)
#'   file. If NULL, if NULL assumes that the plate index is at the end of the
#' .fastq headers the R1 and R3 files.
#' @param R3 character. File name for R3 (read two) file
#' @param indices character. Vector of indices to search for in read files.
#' @param outfile_prefix character, default plate_split. prefix to be appended
#'   to each resulting .fastq file. A vector can be provided if more than one
#'   R1/R2/R3 files are also provided.
#' @param sample_names character. Vector of names for output files corresponding
#'   in order to the provided indices. Prefix and indices will in names will be
#'   replaced with these names. If more than one R1/R2/R3 file is provided, a
#'   name for each index in each plate should be provided as a list in the same
#'   order as the R1/R2/R3 files (for example, \code{list(c("plate1_sample1",
#'   "plate1_sample2"), c("plate2_sample1", "plate2_smaple2"))})
#' @param par numeric, default 1. Number of cores to use for the demultiplexing.
#'   Only used if more than one R1/R2/R3 file are provided.
#'
#' @return Generates files split by index in the directory of the R1 file(s),
#'   named outfile_prefix_RN_INDEX.fastq, where N is 1, 2, 3 three (for the R1,
#'   R2, and R3 files). In R, returns a list containing the paths to each output
#'   file.
#'
#' @export
#'
#' @author William Hemstrom
#' @author Michael Miller
plate_split <- function(R1, R2 = NULL, R3, indices, outfile_prefix = "plate_split", sample_names = NULL,
                        par = 1){

  #============sanity checks========
  msg <- character()

  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }
  
  if(is.null(R2)){
    index_in_header <- TRUE
  }
  else{
    index_in_header <- FALSE
  }
  
  if(!is.character(outfile_prefix)){
    msg <- c(msg, "outfile_prefix must be a character vector of length 1.\n")
  }
  else if(length(outfile_prefix) != 1){
    msg <- c(msg, "outfile_prefix must be a character vector of length 1.\n")
  }
  else if(nchar(outfile_prefix) == 0){
    msg <- c(msg, "outfile_prefix cannot be empty ('').\n")
  }
  
  if(!.paired_length_check(R1, R2) | !.paired_length_check(R1, R3)){
    msg <- c(msg, "R1, R2, and R3 must be of equal length.\n")
  }
  
  if(!is.null(sample_names)){
    if(length(R1) == 1){
      if(!.paired_length_check(indices, unlist(sample_names))){
        msg <- c(msg, "Provided sample_names must be the same length as indices.")
      }
    }
    else{
      if(any(lapply(sample_names, function(x){!.paired_length_check(x, indices)}))){
        msg <- c(msg, "All provided sample_names elements must be the same length as indices.")
      }
      if(!.paired_length_check(sample_names, R1)){
        msg <- c(msg, "sample_names must be equal to the number of R1/R2/R3 files.")
      }
    }
  }
  
  if(length(R1) != 1 & length(outfile_prefix) != 1){
    if(!.paired_length_check(R1, outfile_prefix)){
      msg <- c(msg, "Length of outfile_prefix must either be 1 or equal to the number of R1/R2/R3 files.")
    }
  }

  if(length(msg) > 0){
    stop(msg)
  }
  
  if(par != 1 & length(R1) == 1){
    par <- 1
  }
  
  R1 <- normalizePath(R1)
  if(!is.null(R2)){R2 <- normalizePath(R2)}
  R3 <- normalizePath(R3)
  
  dir <- dirname(R1)

  #============execute=======
  if(length(R1) != 1 & length(outfile_prefix) == 1){
    disambig_prefix <- paste0("file_", 1:length(R1), "_", outfile_prefix)
  }
  else{
    disambig_prefix <- outfile_prefix
  }
  
  disambig_prefix <- paste0(dir, "/", disambig_prefix)
  
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)
  
  # run
  output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    if(index_in_header){
      script <- .fetch_a_script("BarcodeSplitList2Files.pl", "perl")
      cmd <- paste0("perl ", script, " ", R1[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", disambig_prefix[q])
    }
    else{
      script <- .fetch_a_script("BarcodeSplitList3Files.pl", "perl")
      cmd <- paste0("perl ", script, " ", R1[q], " ", R2[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", disambig_prefix[q])
    }
    system(cmd)
    
    
    file_names <- list(R1 = paste0(disambig_prefix[q], "_R1_", indices, ".fastq"),
                       R2 = paste0(disambig_prefix[q], "_R2_", indices, ".fastq"),
                       R3 = paste0(disambig_prefix[q], "_R3_", indices, ".fastq"))
    file_names
  }
  
  file_names <- list(R1 = unlist(purrr::map(output, "R1")),
                     R2 = unlist(purrr::map(output, "R2")),
                     R3 = unlist(purrr::map(output, "R3")))
  
  
  # rename if requested
  if(!is.null(sample_names)){
    sample_names <- unlist(sample_names)
    .rename_files(file_names$R1, paste0(sample_names, "_R1", ".fastq"))
    file_names$R1 <-paste0(sample_names, "_R1", ".fastq")
    
    .rename_files(file_names$R3, paste0(sample_names, "_R3", ".fastq"))
    file_names$R3 <-paste0(sample_names, "_R3", ".fastq")
    
    if(!index_in_header){
      .rename_files(file_names$R2, paste0(sample_names[i], "_R2", ".fastq"))
      file_names$R2 <-paste0(sample_names[i], "_R2", ".fastq")
    }
  }
  
  if(index_in_header){
    file_names$R2 <- NULL
  }
  
  return(file_names)
}

#' Split multiplexed reads with sample barcodes.
#'
#' Splits paired-end sequencing .fastq files containing reads from multiple
#' different plates with unique plate barcodes into different fastq files for
#' each barcode. .fastq files must have barcode identifiers starting at the 5th
#' base of sequence data.
#'
#' Currently, this only supports "strict" barcode matching without *any*
#' sequence mismatches. Supports single- or paired-end sequence data.
#'
#' @param R1 character. File name for R1 (read one) file.
#' @param R2 character, default NULL. File name for R2 (read one) file. If NULL,
#'   assumes data is single-end.
#' @param barcodes character. Vector of barcodes to search for in read files.
#' @param outfile_prefix character. Prefix to be appended to each resulting
#'   .fastq file. A vector can be provided if more than one R1/R2 files are also
#'   provided.
#' @param sample_names character. Vector of names for output files corresponding
#'   in order to the provided barcodes. Prefix and barcodes will in names will
#'   be replaced with these names. If more than one R1/R2 file is provided, a
#'   name for each barcode in each plate should be provided as a list in the
#'   same order as the R1/R2 files (for example, \code{list(c("plate1_sample1",
#'   "plate1_sample2"), c("plate2_sample1", "plate2_smaple2"))})
#' @param stacks_header logical, default TRUE. If TRUE, will fix fastq headers
#'   to be consistent with those expected by stacks (a unique header ending in
#'   either /1 or /2 for read one and two, respectively). If FALSE, headers are
#'   not changed. Only applicable to paired-end sequence data.
#' @param par numeric, default 1. Number of cores to use for the demultiplexing.
#'   Only used if more than one R1/R2 file are provided.
#'
#' @return Generates files split by barcode in the directory of the R1 file(s),
#'   named SAMPLENAME_RN.fastq or outfile_prefix_RN_BARCODE.fastq, if provided
#'   with sample names or not, respectively, where N is A or B (for the reads
#'   with and without barcodes, respectively). In R, returns a list containing
#'   the paths to each output file.
#'   
#' @export
#'
#' @author William Hemstrom
#' @author Michael Miller
demultiplex <- function(R1, R2 = NULL, barcodes, outfile_prefix = "alignR", 
                        sample_names = NULL, stacks_header = TRUE){
  #============sanity checks========
  msg <- character()

  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }
  
  if(!.paired_length_check(R1, R2)){
    msg <- c(msg, "R1 and R2 must be of equal length.\n")
  }
  
  if(!is.null(sample_names)){
    if(length(R1) == 1){
      if(!.paired_length_check(barcodes, unlist(sample_names))){
        msg <- c(msg, "Provided sample_names must be the same length as barcodes.")
      }
    }
    else{
      if(any(lapply(sample_names, function(x){!.paired_length_check(x, barcodes)}))){
        msg <- c(msg, "All provided sample_names elements must be the same length as barcodes.")
      }
      if(!.paired_length_check(sample_names, R1)){
        msg <- c(msg, "sample_names must be equal to the number of R1/R2 files.")
      }
    }
  }
  
  if(length(R1) != 1 & length(outfile_prefix) != 1){
    if(!.paired_length_check(R1, outfile_prefix)){
      msg <- c(msg, "Length of outfile_prefix must either be 1 or equal to the number of R1/R2 files.")
    }
  }
  

  if(length(msg) > 0){
    stop(msg)
  }
  
  if(par != 1 & length(R1) == 1){
    par <- 1
  }
  
  R1 <- normalizePath(R1)
  if(!is.null(R2)){R2 <- normalizePath(R2)}
  
  dir <- dirname(R1)
  #============execute============
  
  if(length(R1) != 1 & length(outfile_prefix) == 1){
    disambig_prefix <- paste0("file_", 1:length(R1), "_", outfile_prefix)
  }
  else{
    disambig_prefix <- outfile_prefix
  }
  
  disambig_prefix <- paste0(dir, "/", disambig_prefix)
  
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)
  
  # run
  output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    
    # run
    if(!is.null(R2)){
      script <- .fetch_a_script("BarcodeSplitListBestRadPairedEnd.pl", "perl")
      cmd <-paste0("perl ", script, " ", R1[q], " ", R2[q], " ", 
                   paste0(barcodes, collapse = ","), " ", disambig_prefix[q],
                   " ", ifelse(stacks_header, 1, 0))
      system(cmd)
    }
    else{
      script <- .fetch_a_script("BarcodeSplitListBestRadSingleEnd.pl", "perl")
      cmd <-paste0("perl ", script, " ", R1[q], " ", paste0(barcodes, collapse = ","), " ", disambig_prefix[q])
      system(cmd)
    }
    
    file_names <- list(RA = paste0(disambig_prefix[q], "_RA_", barcodes, ".fastq"),
                      RB = NULL)
    if(!is.null(R2)){
      file_names$RB <- paste0(disambig_prefix[q], "_RB_", barcodes, ".fastq")
    }
    
    file_names
  }
  
  #============cleanup=============
  file_names <- list(RA = unlist(purrr::map(output, "RA")),
                     RB = unlist(purrr::map(output, "RB")))
  if(is.null(R2)){
    file_names$RB <- NULL
  }

  # rename
  if(!is.null(sample_names)){
    .rename_files(file_names$RA, paste0(sample_names, "_RA", ".fastq"))
    file_names$RA <- paste0(sample_names, "_RA", ".fastq")
    
    if(!is.null(R2)){
      .rename_files(file_names$RB, paste0(sample_names, "_RB", ".fastq"))
      file_names$RB <- paste0(sample_names, "_RB", ".fastq")
    }
  }
  
  if(is.null(file_names[[2]])){file_names[[2]] <- NULL}
  
  return(file_names)
}








