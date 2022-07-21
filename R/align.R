
#' Filter and align single or paired-end reads to a reference genome.
#'
#' Aligns split, single- or paired-end reads to a reference genome using the bwa mem
#' algorithm, them filters the resulting alignments using a few different
#' approaches via samtools.
#' 
#' Paired-end reads are filtered following alignment by first removing PCR
#' duplicates using `samtools fixmate` and `samtools markdup`, then poorly
#' mapped reads using `samtools view` with `-q 5`, and then improperly paired
#' reads using `samtools view` with  `-f 0x2`. The full alignment script can be
#' accessed using `system.file("shell", "run_align.sh", "alignR")`.
#' 
#' Single-end reads are not filtered following alignment, although filtering can
#' (and usually should) be done while calling genotypes or doing other 
#' down-stream analyses.cThe full alignment script can be
#' accessed using `system.file("shell", "run_align_single.sh", "alignR")`.
#' 
#' @param RA_fastqs character. Vector of filepaths for the RA (read one) files.
#' @param RB_fastqs character or NULL, default NULL. Vector of filepaths for the
#'   RB (read one) files.
#' @param reference character. Filepath for the reference genome to use for the
#'   alignment.
#' @param par numeric, default 1. Number of cores to use for the alignments.
#' 
#' @author William Hemstrom
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
  RB_fastqs <- normalizePath(RB_fastqs)
  
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
