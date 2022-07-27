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
#'   .fastq headers the R1 and R3 files.
#' @param R3 character. File name for R3 (read two) file
#' @param indices character. Vector of indices to search for in read files.
#' @param outfile_prefix character, default plate_split. prefix to be appended
#'   to each resulting .fastq file.
#'
#' @export
#'
#' @author William Hemstrom
#' @author Michael Miller
plate_split <- function(R1, R2 = NULL, R3, indices, outfile_prefix = "plate_split"){

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

  if(length(msg) > 0){
    stop(msg)
  }
  
  R1 <- normalizePath(R1)
  if(!is.null(R2)){R2 <- normalizePath(R2)}
  R3 <- normalizePath(R3)

  #============execute=======
  if(index_in_header){
    script <- .fetch_a_script("BarcodeSplitList2Files.pl", "perl")
    cmd <- paste0("perl ", script, " ", R1, " ", R3, " ", paste0(indices, collapse = ","), " ", outfile_prefix)
  }
  else{
    script <- .fetch_a_script("BarcodeSplitList3Files.pl", "perl")
    cmd <- paste0("perl ", script, " ", R1, " ", R2, " ", R3, " ", paste0(indices, collapse = ","), " ", outfile_prefix)
  }
  system(cmd)
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
#' @param outfile_prefix character. prefix to be appended to each resulting
#'   .fastq file.
#' @param sample_names character. Vector of names for output files corresponding
#'   in order to the provided barcodes. Prefix and barcodes will in names
#'   will be replaced with these names.
#' @param stacks_header logical, default TRUE. If TRUE, will fix fastq headers
#'   to be consistent with those expected by stacks (a unique header ending in
#'   either /1 or /2 for read one and two, respectively). If FALSE, headers are
#'   not changed. Only applicable to paired-end sequence data.
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

  if(length(msg) > 0){
    stop(msg)
  }
  
  R1 <- normalizePath(R1)
  if(!is.null(R2)){R2 <- normalizePath(R2)}

  #============execute=======

  # run
  if(!is.null(R2)){
    script <- .fetch_a_script("BarcodeSplitListBestRadPairedEnd.pl", "perl")
    cmd <-paste0("perl ", script, " ", R1, " ", R2, " ", 
                 paste0(barcodes, collapse = ","), " ", outfile_prefix,
                 " ", ifelse(stacks_header, 1, 0))
    system(cmd)
  }
  else{
    script <- .fetch_a_script("BarcodeSplitListBestRadSingleEnd.pl", "perl")
    cmd <-paste0("perl ", script, " ", R1, " ", paste0(barcodes, collapse = ","), " ", outfile_prefix)
    system(cmd)
  }

  # rename
  if(!is.null(sample_names)){
    for(i in 1:length(barcodes)){
      cmdRA <- paste0("mv ", outfile_prefix, "_RA_", barcodes[i], ".fastq", " ", sample_names[i], "_RA", ".fastq")
      system(cmdRA)
      
      if(!is.null(R2)){
        cmdRB <- paste0("mv ", outfile_prefix, "_RB_", barcodes[i], ".fastq", " ", sample_names[i], "_RB", ".fastq")
        system(cmdRB)
      }
    }
  }

}








