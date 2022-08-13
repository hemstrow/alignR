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
demultiplex <- function(R1, R2 = NULL, R3 = NULL, 
                        barcodes = NULL, indices = NULL, 
                        barcode_locations, 
                        outfile_prefix = "alignR", 
                        sample_names = NULL, 
                        stacks_header = TRUE, 
                        par = 1){
  #============sanity checks========
  msg <- character()

  # installed
  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }
  
  # R1/R2/R3 lengths
  if(!.paired_length_check(R1, R2)){
    msg <- c(msg, "R1, R2, and R3 (if provided) must be of equal length.\n")
  }
  else if(!.paired_length_check(R1, R3)){
    msg <- c(msg, "R1, R2 (if provided), and R3 must be of equal length.\n")
  }
  
  
  # sample names for the list of DFs
  .check_sample_names <- function(sample_names, barcodes, indices){
    msg <- character()
    
    # sample names and prefixes
    if(!is.null(sample_names)){
      
      # one or the other
      if(xor(is.null(barcodes), is.null(indices))){
        
        # correct number of rows?
        if(nrow(sample_names) != length(barcodes) + length(indices)){
          msg <- c(msg, "Number of rows in sample_names must equal the number of provided barcodes/indices.\n")
        }
        
        # have 2 columns?
        if(ncol(sample_names) != 2){
          msg <- c(msg, "If only one of barcodes or indices are provided, sample_names must have two columns.\n")
        }
        # have all the barcodes acounted for?
        else{
          if(!is.null(barcodes)){
            bad.barcodes <- which(!barcodes %in% sample_names[,2])
            if(length(bad.barcodes) != 0){
              msg <- c(msg, paste0("Some barcodes not found in sample_names data.frame: \n\t", paste0(barcodes[bad.barcodes], collapse = "\n\t"), "\n"))
            }
          }
          else{
            bad.indices <- which(!indices %in% sample_names[,2])
            if(length(bad.indices) != 0){
              msg <- c(msg, paste0("Some indices not found in sample_names data.frame: \n\t", paste0(indices[bad.indices], collapse = "\n\t"), "\n"))
            }
          }
          
        }
      }
      if(!is.null(barcodes) & !is.null(indices)){
        
        # have correct number of rows?
        if(nrow(sample_names) != length(barcodes) * length(indices)){
          msg <- c(msg, "Number of rows in sample_names must equal the number of possible combinations for the provided barcodes/indices.\n")
        }
        
        # have the correct number of columns?
        if(ncol(sample_names) != 3){
          msg <- c(msg, "Three columns (names, barcodes, indices, in that order) are needed in sample_names.\n")
        }
        # check all barcodes and index combos are accounted for
        else{
          needed_rows <- expand.grid(barcodes, indices)
          needed_rows <- do.call(paste, c(needed_rows[,1:2, drop = FALSE], sep = "\t"))
          have_rows <- do.call(paste, c(sample_names[,2:3, drop = FALSE], sep = "\t"))
          missing_combos <-  which(!needed_rows %in% have_rows)
          if(length(missing_combos) != 0){
            msg <- c(msg, paste0("Some barcode/index combinations not found in sample_names:\n\t", paste0(needed_rows[missing_combos], collapse = "\n\t"), "\n"))
          }
        }
      }
    }
    
    return(msg)
  }
  
  if(!is.data.frame(sample_names) & is.list(sample_names)){
    sn_checks <- vector("list", length(sample_names))
    for(i in 1:length(sample_names)){
      sn_checks[[i]] <- .check_sample_names(x, barcodes, indices)
      if(length(sn_checks[[i]]) != 0){
        sn_checks[[i]] <- paste0("Problem with sample_names df ", i, " : ", sn_checks[[i]])
      }
    }
    sn_checks <- unlist(sn_checks)
    msg <- c(msg, sn_checks)
  }
  
  
  
  # barcode options
  good.options <- c("header", "read_start", "R2")
  bad_barcode_locations <- which(!barcode_locations %in% good.options)
  if(lenth(bad_barcode_locations) != 0){
    msg <- c(msg, paste0("Unaccepted barcode locations: ", paste0(barcode_locations[bad_barcode_locations], collapse = "  "), "\n"))
  }
  if(all(c("header", "read_start") %in% barcode_locations)){
    msg <- c(msg, "Cannot demultiplex reads with barcodes both in the header and at the start of reads.")
  }
  
  # check that we have the right barcodes/indices and reads
  if("header" %in% barcode_locations | "read_start" %in% barcode_locations){
    if(is.null(barcodes)){
      msg <- c(msg, "barcodes must be provided if barcodes are located in the header or at the start of reads.\n")
    }
  }
  if("R2" %in% barcode_locations){
    if(is.null(indices)){
      msg <- c(msg, "indices must be provided if barcodes are located in the R2 reads.\n")
    }
    if(is.null(R2)){
      msg <- c(msg, "Path(s) to R2 reads must be provided if barcodes are located in R2 reads.\n")
    }
  }
  
  # stop if bad
  if(length(msg) > 0){
    stop(msg)
  }
  
  # check par
  if(!(!is.null(indices) & !is.null(barcodes))){
    if(length(indices) + length(barcodes) != 1){
      par <- 1
    }
  }
  else if(length(indices) == 1){
    par <- 1
  }
  
  R1 <- normalizePath(R1)
  if(!is.null(R2)){R2 <- normalizePath(R2)}
  if(!is.null(R3)){R3 <- normalizePath(R3)}
  
  dir <- dirname(R1)
  has_prefix <- FALSE
  #============split by R2========
  
  if("R2" %in% barcode_locations){
    
    # prep
    ps_par <- min(par, length(R1))
    
    cl <- parallel::makePSOCKcluster(ps_par)
    doParallel::registerDoParallel(cl)
    
    if(is.null(R2)){
      script <- .fetch_a_script("demulti_single_end_R2_barcodes.pl", "perl")
    }
    else{
      script <- .fetch_a_script("demulti_paired_end_R2_barcodes.pl", "perl")
    }
    
    
    R1_base <- tools::file_path_sans_ext(basename(R1))
    R1_base <- gsub("R1", "", R1_base)
    R1_base <- gsub("__", "_", R1_base)
    
    # run
    output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {
      if(is.null(R2)){
        cmd <- paste0("perl ", script, " ", R1[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", paste0(outfile_prefix, "_", R1_base[q]))
      }
      else{
        cmd <- paste0("perl ", script, " ", R1[q], " ", R2[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", paste0(outfile_prefix, "_", R1_base[q]))
      }
    }
    
    # wrap
    parallel::stopCluster(cl)
    
    if(any(c("headers", "read_start") %in% barcode_locations)){
      outputs <- expand.grid(R1_base, indices)
      R1 <- paste0(outfile_prefix, "_", outputs[,1], "_R1_", outputs[,2], ".fastq")
      R2 <- NULL
      R3 <- NULL
      has_prefix <- TRUE
    }
  }
  
  #============split by header/start of read============
  
  if(any(c("header", "read_start") %in% barcode_locations)){
    
    par <- min(par, length(R1))
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)
    
    R1_base <- tools::file_path_sans_ext(basename(R1))
    R1_base <- gsub("R1", "", R1_base)
    R1_base <- gsub("__", "_", R1_base)
    
    output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {
      
      # run
      if(!is.null(R3)){
        script <- .fetch_a_script("demulti_single_end_single_file.pl", "perl")
        cmd <-paste0("perl ", 
                     script, " ", 
                     R1[q], " ",
                     paste0(barcodes, collapse = ","), " ", 
                     paste0(ifelse(has_prefix, "", paste0(outfile_prefix), "_"), R1_base[q])," ", 
                     ifelse(stacks_header, 1, 0),
                     ifelse("header" %in% barcode_locations, " 1", " 0"))
      }
      else{
        script <- .fetch_a_script("demulti_paired_end_single_file.pl", "perl")
        cmd <-paste0("perl ", 
                     script, " ", 
                     R1[q], " ", 
                     R3[q], " ", 
                     paste0(barcodes, collapse = ","), " ", 
                     paste0(ifelse(has_prefix, "", paste0(outfile_prefix), "_"), R1_base[q])," ", 
                     ifelse(stacks_header, 1, 0),
                     ifelse("header" %in% barcode_locations, " 1", " 0"))
        
      }
      
      system(cmd)
      
      parallel::stopCluster(cl)
    }
  }

  #============rename and finish=============


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








