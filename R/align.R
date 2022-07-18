
#' Filter and align paired-end reads to a reference genome.
#'
#'
align_to_reference <- function(RA_fastqs, RB_fastqs, reference, par = 1){
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
  #==============prepare to run=================
  script <- .fetch_a_script("run_align.sh", "shell")

  # prepare file handles
  fastqs <- data.frame(RA = RA_fastqs, RB = RB_fastqs)
  fastqs$RA <- gsub("\\.fastq$", "", fastqs$RA)
  fastqs$RB <- gsub("\\.fastq$", "", fastqs$RB)

  # register parallel architecture
  cl <- snow::makeSOCKcluster(par)
  doSNOW::registerDoSNOW(cl)

  # divide up into ncore chunks
  iters <- length(RA_fastqs)
  it_par <- (1:iters)%%par
  chunks <- .smart_split(fastqs, it_par)

  # prepare reporting function
  progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
  opts <- list(progress=progress)


  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .options.snow = opts, .packages = c("alignR")
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      cmd <- paste0("bash ", script, " ", fastqs$RA[i], " ", fastqs$RB[i], " ", reference)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)
  doSNOW::registerDoSNOW()
  gc();gc()
}
