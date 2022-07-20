#' Extract SNP genotypes from aligned bamfiles using ANGSD.
#'
#' Uses the ANGSD (Analysis of Next Gen Sequencing Data) software to genotype
#' individuals across a set of bamfiles using several different filtering
#' options or genotyping methods. This function is a \emph{ wrapper for ANGSD}.
#'
#' @param bamfiles character. Vector of filepaths to the bamfiles containing
#'   data for the individuals to genotype.
#' @param outfile character. Prefix for the outfiles. Can include a filepath.
#' @param minInd numeric, default \code{floor(length(bamfiles)/2)}. Minimum number of
#'   individuals a locus must be sequenced in in order to call genotypes.
#' @param genotyper character, default "SAMtools". Name of the genotyper to use.
#'   Options: \itemize{\item{SAMtools} \item{GATK} \item{SOAPsnp} \item{phys}
#'   \item{sample: } randomly sample reads to make the genotypes (therefore
#'   posteriors are either 0, .5, or 1).} See ANGSD documentation for details.
#' @param SNP_pval numeric, default 0.00000001. Minimum p-value that must be met
#'   \emph{for the probability that a given locus is a SNP} for a locus to be
#'   genotyped.
#' @param doGeno character, default "NN". Vector of styles for which genotypes
#'   will be printed. Options: \itemize{\item{major_minor: } Print the major and
#'   minor alleles. \item{numeric: } Print alleles as 0, 1, 2, or -1 for the
#'   homozygous major, hetoerzygote, homozygous minor, or missing data,
#'   respectively. \item{NN: } Print alleles as paired nucleotides (AA, AC, CC,
#'   for example). NN denotes missing data. \item{all_posteriors: } Prints
#'   posterior probabilities for \emph{all} genotypes at each locus.
#'   \item{called_posterior: } Prints the posterior probability for \emph{only
#'   the called genotype}. \item{SNP_genotype_posteriors: } Prints the posterior
#'   probabilities for \emph{the three genotypes possible for the major and
#'   minor allele} at each locus only.}
#' @param postCutoff numeric, default 0.95. The minimum posterior probability
#'   required for the most likely genotype for the genotype to be printed at for
#'   a given individual at a given locus.
#' @param minQ numeric, default 20. The minimum Phred \emph{sequencing} quality
#'   score needed for a given sequenced base on a single read to be considered
#'   during genotyping. The default, 20, corresponds to 99\% accuracy.
#' @param minMapQ numeric, default 20. The minimum Phred \emph{mapping} quality
#'   score needed for a given sequenced base on a single read to be considered
#'   during genotyping. The default, 20, corresponds to 99\% accuracy.
#' @param par numeric, default 1. Number of cores to allow ANGSD to use for
#'   genotyping.
#'
#' @author William Hemstrom
#' @author Michael Miller
#'
#' @export
genotype_bams <- function(bamfiles,
                          outfile,
                          minInd = floor(length(bamfiles)/2),
                          genotyper = "SAMtools",
                          SNP_pval = 0.00000001,
                          doGeno = "NN",
                          postCutoff = 0.95,
                          minQ = 20,
                          minMapQ = 20,
                          par = 1){
  #=============sanity checks===================
  msg <- character()

  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  if(!.check_system_install("angsd")){
    msg <- c(msg, "No angsd install located on system path.\n")
  }

  if(length(msg) > 0){
    stop(msg)
  }

  # check that the bams exist and are indexed.
  if(all(file.exists(bamfiles))){
    index_files <- paste0(bamfiles, ".bai")
    indexed <- file.exists(index_files)
    if(any(!indexed)){
      cat(sum(!indexed), "bam files not indexed. Indexing with 'samtools faidx'.\n")

      if(!.check_system_install("samtools")){
        stop("No samtools install located on system path.\n")
      }

      for(i in 1:sum(!indexed)){
        system(paste0("samtools faidx ", bamfiles[!indexed][i]))
      }
    }
  }
  else{
    msg <- c(msg, paste0("Some bamfiles not located: ",
                         paste0(bamfiles[!file.exists(bamfiles)], collapse = ", "),
                         ".\n"))
  }

  if(is.character(read)){
    missing_doGeno <- !read %in% doGeno
    if(any(missing_doGeno)){
      warning("Some read returns requested but not included in doGeno. Adding these to doGeno")
      doGeno <- c(doGeno, read[which(missing_doGeno)])
    }
  }

  # figure out genotyper code
  genotyper_table <- data.frame(matrix(c(1, "SAMtools",
                                         2, "GATK",
                                         3, "SOAPsnp",
                                         4, "SYK",
                                         5, "phys",
                                         6, "sample"),
                                       ncol = 2, byrow = TRUE))
  if(!genotyper %in% genotyper_table[,2]){
    msg <- c(msg, paste0("Genotyper ", genotyper, " not available. Is this a typo?\n"))
  }
  else{
    genotyper <- genotyper_table[match(genotyper, genotyper_table[,2]),1]
  }

  # figure out doGeno
  doGeno_table <- data.frame(number = c(1, 2, 4, 8, 16, 32),
                             word = c("major_minor", "numeric", "NN", "all_posteriors",
                                      "called_posterior", "SNP_genotype_posteriors"))
  bad_doGeno <- !doGeno %in% doGeno_table[,2]
  if(any(bad_doGeno)){
    msg <- c(msg, paste0("doGeno option(s) ", paste0(doGeno[bad_doGeno], collapse = ", "),
                         " not available. Is this a typo?\n"))
  }
  else{
    doGeno <- doGeno_table[match(doGeno, doGeno_table[,2]),1]
    doGeno <- sum(doGeno)
  }

  if(length(msg) > 0){
    stop(msg)
  }

  #==============prepare to run=================
  browser()
  old.scipen <- options("scipen")
  options(scipen = 999)
  script <- .fetch_a_script("angsd_genotypes.sh", "shell")

  # save the bamlist
  write(bamfiles, paste0(outfile, "_bamlist.txt"), ncolumns = 1)

  # compose command and run
  cmd <- paste0("bash ", script, " ",
                paste0(c(
                  paste0(outfile, "_bamlist.txt"),
                  minInd,
                  genotyper,
                  SNP_pval,
                  doGeno,
                  postCutoff,
                  minQ,
                  minMapQ,
                  outfile,
                  par), collapse = " ")
                )

  system(cmd)

  options(scipen = old.scipen)
}
