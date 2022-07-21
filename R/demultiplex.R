#' Split paired-end .fastq files with reads from different plates.
#'
#' Splits paired-end sequencing .fastq files containing reads from multiple
#' different plates with unique plate barcodes into different fastq files for
#' each plate barcode. 
#' 
#' The R1 and R3 .fastq files should contain the actual forward and reverse
#' reads. The R2 .fastq file should contain only the plate barcodes. As is
#' standard, the reads are assumed to be in the same order for each of the
#' .fastq files.
#'
#' @param R1 character. File name for R1 (read one) file.
#' @param R2 character. File name for R2 (barcode read) file.
#' @param R3 character. File name for R3 (read two) file
#' @param barcodes character. Vector of barcodes to search for in read files.
#' @param outfile_prefix character. prefix to be appended to each resulting
#'   .fastq file.
#'
#' @export
#'
#' @author William Hemstrom
#' @author Michael Miller
plate_split <- function(R1, R2, R3, barcodes, outfile_prefix = ""){

  #============sanity checks========
  msg <- character()

  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }

  if(length(msg) > 0){
    stop(msg)
  }
  
  R1 <- normalizePath(R1)
  R2 <- normalizePath(R2)
  R3 <- normalizePath(R3)

  #============execute=======
  script <- .fetch_a_script("BarcodeSplitList3Files.pl", "perl")
  cmd <- paste0("perl ", script, " ", R1, " ", R2, " ", R3, " ", paste0(barcodes, collapse = ","), " ", outfile_prefix)
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
#' sequence mismatches.
#'
#' @param R1 character. File name for R1 (read one) file.
#' @param R2 character. File name for R2 (read one) file.
#' @param barcodes character. Vector of barcodes to search for in read files.
#' @param outfile_prefix character. prefix to be appeneded to each resulting
#'   .fastq file.
#' @param sample_names character. Vector of names for output files corresponding
#'   in order to the provided barcodes. Prefix and barcodes will in names
#'   will be replaced with these names.
#'
#' @export
#'
#' @author William Hemstrom
#' @author Michael Miller
demultiplex <- function(R1, R2, barcodes, outfile_prefix = "alignR", sample_names = NULL){
  #============sanity checks========
  msg <- character()

  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }

  if(length(msg) > 0){
    stop(msg)
  }
  
  R1 <- normalizePath(R1)
  R2 <- normalizePath(R2)

  #============execute=======

  # run
  script <- .fetch_a_script("BarcodeSplitListBestRadPairedEnd.pl", "perl")
  cmd <-paste0("perl ", script, " ", R1, " ", R2, " ", paste0(barcodes, collapse = ","), " ", outfile_prefix)
  system(cmd)

  # rename
  if(!is.null(sample_names)){
    for(i in 1:length(barcodes)){
      cmdRA <- paste0("mv ", outfile_prefix, "_RA_", barcodes[i], ".fastq", " ", sample_names[i], "_RA", ".fastq")
      cmdRB <- paste0("mv ", outfile_prefix, "_RB_", barcodes[i], ".fastq", " ", sample_names[i], "_RB", ".fastq")
      system(cmdRA)
      system(cmdRB)
    }
  }

}








