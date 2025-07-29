#' Split multiplexed reads via barcodes.
#'
#' Splits paired-end or single-end, single or dual-indexed sequencing fastq
#' files containing reads from multiple individuals.
#'
#' Currently, this only supports "strict" barcode matching without \emph{any}
#' sequence mismatches. Supports single- or paired-end sequence data that is
#' either dual- or single-indexed. Barcodes can be located either in a separate
#' sequencing file (usually identified by the string "R2" somewhere in the file
#' names), in the fastq headers (at the end of the header, following the final
#' ":"), or at the start of either the forward or reverse read.
#'
#' @param R1 character vector. File name(s) for R1 (forward read) file(s).
#' @param R2 character vector, default NULL. File name(s) for R2 (barcode) read
#'   file(s). Needed if barcodes are located on a separate read.
#' @param R3 character vector, default NULL. File name(s) for the R3 (reverse
#'   read) file(s). If NULL, assumes the data is single-end.
#' @param barcodes character. Vector of \emph{header or start-of-read} barcodes
#'   or file path to file containing such with one barcode on each line. Needed
#'   if either the \code{read_start} or \code{header} \code{barcode_location}
#'   options are selected.
#' @param indices character. Vector of \emph{R2, seperate read} barcodes or file
#'   path to file containing such with one barcode on each line. Needed if
#'   the \code{R2} \code{barcode_location} option is selected.
#' @param barcode_locations character. String indicating the locations of
#'   barcodes. Options: \itemize{\item{\code{read_start}: } barcodes located at
#'   the start of each read. \item{\code{header}: } barcodes located in the
#'   fastq header for each read, after the final ":". \item{\code{R2}: }
#'   barcodes located in separate fastq files (usually named ending in "R2").
#'   Reads in these files contain \emph{only} the barcodes.} \code{R2} can be
#'   selected alongside either \code{header} or \code{read_start}. The latter
#'   two cannot currently be selected together.
#' @param outfile_prefix character. Prefix to be appended to each resulting
#'   .fastq file.
#' @param sample_names Optional data.frame or list of data.frames. A data.frame
#'   for each input R1/R2/R3 file (so if 3 of each are provided, this should be
#'   a list of 3 data frames) with either two or three columns containing:
#'   \enumerate{\item Unique sample IDs \item barcodes provided to
#'   \code{indices} if \code{R2} \code{barcode_location} is selected, barcodes
#'   provided to \code{barcodes} otherwise. \item barcodes provided to
#'   \code{barcodes} if both \code{R2} and either \code{read_start} or
#'   \code{header} were selected for \code{barcode_location}.} In essence, the
#'   second column will hold indices if provided, which bumps barcodes into
#'   column three. The second column will otherwise hold barcodes. If multiple
#'   R1/R2/R3 files are provided, \emph{data.frames will be evaluated and
#'   matched to samples in order} (the first R1/R2/R3 file will correspond to
#'   the first data.frame and so on). \emph{Each barcode, index, or
#'   barcode/index combination must be present in each data.frame!}
#' @param stacks_header logical, default TRUE. If TRUE, will fix fastq headers
#'   to be consistent with those expected by stacks (a unique header ending in
#'   either /1 or /2 for read one and two, respectively). If FALSE, headers are
#'   not changed. Only applicable to paired-end sequence data.
#' @param par numeric, default 1. Number of cores to use for the demultiplexing.
#'   Only used if more than one R1/R2 file are or data is dual-indexed.
#'
#' @return Generates files split by barcodes in the directory of the R1 file(s),
#'   named SAMPLENAME_RN.fastq if \code{sample_names} are provided,
#'   outfile_prefix_RN_BARCODE.fastq, outfile_prefix_RN_INDEX.fastq, or
#'   outfile_prefix_RN_INDEX_RN_BARCODE.fastq if not, where N is either (1/2/3)
#'   or (A/B) for reads split by R2 barcodes or start-of-read/headers,
#'   respectively. R1 reads are forward reads, R2 reads are barcode reads, and
#'   R3 reverse reads. RA reads had barcodes (that were trimmed by
#'   \code{demultiplex}), RB reads did not. Returns a list containing the paths
#'   to each output file.
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

  # check that we have barcodes and they are OK
  if(length(barcodes) == 1){
    if(file.exists(normalizePath(barcodes))){
      barcodes <- readLines(normalizePath(barcodes))
      if(length(barcodes) == 1){
        msg <- c(msg, "Barcodes must either be directly provided or given in a line-seperated file.\n")
      }
    }
  }
  if(!is.null(barcodes)){
    if(any(!unlist(strsplit(barcodes, "")) %in% c("A", "T", "G", "C"))){
      msg <- c(msg, "Unaccepted characters found in barcodes. Accepted characters: A, T, C, G.\n")
    }
  }


  # check that we have indices and they are ok
  if(length(indices) == 1){
    if(file.exists(normalizePath(indices))){
      indices <- readLines(normalizePath(indices))
      if(length(indices) == 1){
        msg <- c(msg, "indices must either be directly provided or given in a line-seperated file.\n")
      }
    }
  }
  if(!is.null(indices)){
    if(any(!unlist(strsplit(indices, "")) %in% c("A", "T", "G", "C"))){
      msg <- c(msg, "Unaccepted characters found in indices. Accepted characters: A, T, C, G.\n")
    }
  }


  # function to check sample names for the list of DFs
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
          needed_rows <- expand.grid(indices, barcodes)
          needed_rows <- paste0(needed_rows[,1], "\t", needed_rows[,2])
          have_rows <- paste0(sample_names[,2], "\t", sample_names[,3])
          missing_combos <-  which(!needed_rows %in% have_rows)
          if(length(missing_combos) != 0){
            msg <- c(msg, paste0("Some barcode/index combinations not found in sample_names:\n\t", paste0(needed_rows[missing_combos], collapse = "\n\t"), "\n"))
          }
        }
      }
    }

    return(msg)
  }

  # issues with sample names?
  if(!is.data.frame(sample_names) & is.list(sample_names)){
    ln_sample_names <- length(sample_names)
    sn_checks <- vector("list", ln_sample_names)
    for(i in 1:ln_sample_names){
      if(!is.data.frame(sample_names[[i]])){stop(c(msg, "sample_names must be a data.frame or list of data.frames.\n"))}
      sn_checks[[i]] <- .check_sample_names(sample_names[[i]], barcodes, indices)
      if(length(sn_checks[[i]]) != 0){
        sn_checks[[i]] <- paste0("Problem with sample_names df ", i, " : ", sn_checks[[i]])
      }
    }
    sn_checks <- unlist(sn_checks)


    if(ln_sample_names != length(R1)){
      sn_checks <- c(sn_checks, "The length of the sample_names list must be equal to the number of R1/R2/R3 files (each data.frame in the list corresponds in order to one R1/R2/R3 file).\n")
    }
    msg <- c(msg, sn_checks)
  }
  else if(!is.null(sample_names)){
    if(!is.data.frame(sample_names)){stop(c(msg, "sample_names must be a data.frame or list of data.frames.\n"))}
    ln_sample_names <- 1
    msg <- c(msg, .check_sample_names(sample_names, barcodes, indices))

    sample_names <- list(sample_names)
  }

  # all names unique?
  if(!is.null(sample_names)){
    all_names <- unlist(purrr::map(sample_names, 1))
    dup_names <- duplicated(all_names)

    if(sum(dup_names) != 0){
      msg <- c(msg, paste0("Some duplicated sample names detected: \n\t", paste0(all_names[dup_names], collapse = "\n\t"), "\n"))
    }
  }



  # barcode options
  good.options <- c("header", "read_start", "R2")
  bad_barcode_locations <- which(!barcode_locations %in% good.options)
  if(length(bad_barcode_locations) != 0){
    msg <- c(msg, paste0("Unaccepted barcode locations: ", paste0(barcode_locations[bad_barcode_locations], collapse = "  "), "\n"))
  }
  if(all(c("header", "read_start") %in% barcode_locations)){
    msg <- c(msg, "Cannot demultiplex reads with barcodes both in the header and at the start of reads.\n")
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

  dirs <- dirname(c(R1, R2, R3))
  if(length(unique(dirs)) != 1){
    msg <- c(msg, "All input fastq files must be in the same directory.\n")
  }

  # stop if bad
  if(length(msg) > 0){
    stop(msg)
  }

  # check par

  R1 <- normalizePath(R1)
  original_R1 <- R1
  if(!is.null(R2)){R2 <- normalizePath(R2)}
  if(!is.null(R3)){R3 <- normalizePath(R3)}

  dir <- dirname(R1[1])
  has_prefix <- FALSE

  if(outfile_prefix != ""){
    outfile_prefix <- paste0(outfile_prefix, "_")
  }
  #============figure out output file names beforehand============

  bns <- basename(R1)
  outfiles <- data.frame(initial_bn = bns)


  # if we have R2 barcodes
  if("R2" %in% barcode_locations){

    outfiles <- expand.grid(outfiles$initial_bn, indices, stringsAsFactors = F) # possible options


    # R1 and R2 (and R3 if paried) files
    outfiles$R1 <- file.path(dir, paste0(outfile_prefix,
                                         tools::file_path_sans_ext(outfiles$Var1), "_R1_",
                                         outfiles$Var2,
                                         ".fastq"))
    outfiles$R2 <- file.path(dir, paste0(outfile_prefix,
                                         tools::file_path_sans_ext(outfiles$Var1), "_R2_",
                                         outfiles$Var2,
                                         ".fastq"))
    if(!is.null(R3)){
      outfiles$R3 <- file.path(dir, paste0(outfile_prefix,
                                           tools::file_path_sans_ext(outfiles$Var1), "_R3_",
                                           outfiles$Var2,
                                           ".fastq"))
    }

    colnames(outfiles)[1:2] <- c("original_R1", "index")

    # if we also have read start/header barcodes, return RA (and RB if paired)
    if(length(barcode_locations) > 1){
      expanded_outfiles <- expand.grid(outfiles$R1, barcodes) # possible options
      expanded_outfiles$original_R1 <- outfiles$original_R1[match(expanded_outfiles$Var1, outfiles$R1)]
      expanded_outfiles$index <- outfiles$index[match(expanded_outfiles$Var1, outfiles$R1)]

      expanded_outfiles$RA <- paste0(tools::file_path_sans_ext(expanded_outfiles$Var1), "_RA_",
                                     expanded_outfiles$Var2, ".fastq")
      if(!is.null(R3)){
        expanded_outfiles$RB <- paste0(tools::file_path_sans_ext(expanded_outfiles$Var1), "_RB_",
                                       expanded_outfiles$Var2, ".fastq")
      }
      colnames(expanded_outfiles)[1:2] <- c("input_R1", "barcode")
    }
  }
  # only read start/header barcodes
  else{
    outfiles <- expand.grid(outfiles$initial_bn, barcodes, stringsAsFactors = F) # possible options
    outfiles$RA <- file.path(dir, paste0(outfile_prefix,
                                         tools::file_path_sans_ext(outfiles$Var1), "_RA_",
                                         outfiles[,2], ".fastq"))
    if(!is.null(R3)){
      outfiles$RB <- file.path(dir, paste0(outfile_prefix,
                                           tools::file_path_sans_ext(outfiles$Var1), "_RB_",
                                           outfiles[,2], ".fastq"))
    }
    colnames(outfiles)[1:2] <- c("original_R1", "barcode")

  }

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

    # run
    output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {
      out <- paste0(file.path(dir, outfile_prefix), basename(tools::file_path_sans_ext(R1)[q]))

      if(is.null(R2)){
        cmd <- paste0("perl ", script, " ", R1[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", out)
      }
      else{
        cmd <- paste0("perl ", script, " ", R1[q], " ", R2[q], " ", R3[q], " ", paste0(indices, collapse = ","), " ", out)
      }

      system(cmd)

      out
    }

    # wrap
    parallel::stopCluster(cl)

    if(length(barcode_locations) != 1){
      R1 <- outfiles$R1
      if(!is.null(R3)){R3 <- outfiles$R3}

      has_prefix <- TRUE
    }
  }

  #============split by header/start of read============

  if(any(c("header", "read_start") %in% barcode_locations)){

    par <- min(par, length(R1))
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)

    if(is.null(R3)){
      script <- .fetch_a_script("demulti_single_end_single_file.pl", "perl")
    }
    else{
      script <- .fetch_a_script("demulti_paired_end_single_file.pl", "perl")
    }


    output <- foreach::foreach(q = 1:length(R1), .inorder = TRUE, .errorhandling = "pass",
                               .packages = "alignR"
    ) %dopar% {

      out <- file.path(dir, paste0(ifelse(has_prefix, "", outfile_prefix), basename(tools::file_path_sans_ext(R1)[q])))

      # run
      if(is.null(R3)){
        cmd <-paste0("perl ",
                     script, " ",
                     R1[q], " ",
                     paste0(barcodes, collapse = ","), " ",
                     out, " ",
                     ifelse(stacks_header, 1, 0),
                     ifelse("header" %in% barcode_locations, " 1", " 0"))
      }
      else{
        cmd <-paste0("perl ",
                     script, " ",
                     R1[q], " ",
                     R3[q], " ",
                     paste0(barcodes, collapse = ","), " ",
                     out, " ",
                     ifelse(stacks_header, 1, 0),
                     ifelse("header" %in% barcode_locations, " 1", " 0"))

      }
      system(cmd)
      out
    }

    parallel::stopCluster(cl)

    if("R2" %in% barcode_locations){
      outfiles <- expanded_outfiles
    }
  }

  #============rename and finish=============

  # rename according to key
  if(!is.null(sample_names)){

    # R2 barcodes
    if("R2" %in% barcode_locations){

      tfn <- vector("list", length = ln_sample_names)
      for(i in 1:ln_sample_names){

        these_files <- outfiles[which(outfiles$original_R1 == basename(original_R1[i])),, drop = F]
        these_names <- sample_names[[i]]

        # in read as well?
        if(length(barcode_locations) > 1){
          colnames(these_names) <- c("ID", "index", "barcode")
          rename_key <- merge(these_names, these_files, by = c("index", "barcode"))
          rename_key$new_RA <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_RA.fastq"))

          file.rename(rename_key$RA, rename_key$new_RA)

          tfn[[i]] <- list(RA = rename_key$new_RA)

          if("RB" %in% colnames(rename_key)){
            rename_key$new_RB <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_RB.fastq"))

            file.rename(rename_key$RB, rename_key$new_RB)
            tfn[[i]]$RB <- rename_key$new_RB
          }
        }

        # just R2?
        else{
          colnames(these_names) <- c("ID", "index")
          colnames(these_files)[1:2] <- c("file_name", "index")
          rename_key <- merge(these_names, these_files, by = "index")
          rename_key$new_R1 <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_R1.fastq"))
          rename_key$new_R2 <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_R2.fastq"))

          file.rename(rename_key$R1, rename_key$new_R1)
          file.rename(rename_key$R2, rename_key$new_R2)

          tfn[[i]] <- list(R1 = rename_key$new_R1, R2 = rename_key$new_R2)

          if("R3" %in% names(rename_key)){
            rename_key$new_R3 <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_R3.fastq"))

            file.rename(rename_key$R3, rename_key$new_R3)
            tfn[[i]]$R3 <- rename_key$new_R3
          }
        }
      }

      file_names <- dplyr::bind_rows(tfn)
    }
    # other barcodes
    else{

      tfn <- vector("list", length = ln_sample_names)
      for(i in 1:ln_sample_names){
        these_files <- outfiles[which(outfiles$original_R1 == basename(original_R1[i])),, drop = F]
        these_names <- sample_names[[i]]

        colnames(these_names) <- c("ID", "barcode")
        colnames(these_files)[1:2] <- c("file_name", "barcode")
        rename_key <- merge(these_names, these_files, by = "barcode")
        rename_key$new_RA <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_RA.fastq"))

        file.rename(rename_key$RA, rename_key$new_RA)

        tfn[[i]] <- list(RA = rename_key$new_RA)

        if("RB" %in% names(rename_key)){
          rename_key$new_RB <- file.path(dir, paste0(outfile_prefix, rename_key$ID, "_RB.fastq"))

          file.rename(rename_key$RB, rename_key$new_RB)
          tfn[[i]]$RB <- rename_key$new_RB
        }
      }

      file_names <- dplyr::bind_rows(tfn)
    }


    file_names <- as.list(file_names)
  }


  # if no key, fix any R1_R1 etc in names
  else{
    if("R2" %in% barcode_locations){

      # R2 + other
      if(length(barcode_locations) > 1){
        R1_R1s <- grep("R1_R1", outfiles$RA)

        if(length(R1_R1s) != 0){
          new_RA <- gsub("R1_R1", "R1", outfiles$RA[R1_R1s])
          file.rename(outfiles$RA[R1_R1s], new_RA)
          outfiles$RA[R1_R1s] <- new_RA

          if(!is.null(R3)){
            new_RB <- gsub("R1_R1", "R1", outfiles$RB[R1_R1s])
            file.rename(outfiles$RB[R1_R1s], new_RB)
            outfiles$RB[R1_R1s] <- new_RB
          }
        }

        file_names <- list(RA = outfiles$RA)
        if("RB" %in% colnames(outfiles)){file_names$RB <- outfiles$RB}
      }

      # RA only
      else{
        R1_R1s <- grep("R1_R1", basename(outfiles$R1))
        if(length(R1_R1s != 0)){
          new_R1 <- file.path(dir, gsub("R1_R1", "R1", basename(outfiles$R1[R1_R1s])))
          file.rename(outfiles$R1[R1_R1s], new_R1)

          new_R2 <- file.path(dir, gsub("R1_R2", "R2", basename(outfiles$R2[R1_R1s])))
          file.rename(outfiles$R2[R1_R1s], new_R2)

          if(!is.null(R3)){
            new_R3 <- file.path(dir, gsub("R1_R3", "R3", basename(outfiles$R3[R1_R1s])))
            file.rename(outfiles$R3[R1_R1s], new_R3)
          }
        }

        file_names <- list(R1 = outfiles$R1, R2 = outfiles$R2)
        if("R3" %in% colnames(outfiles)){
          file_names$R3 <- outfiles$R3
        }
      }
    }
    else{
      file_names <- list(RA = outfiles$RA)
      if("RB" %in% colnames(outfiles)){
        file_names$RB <- outfiles$RB
      }
    }
  }

  return(file_names)
}








