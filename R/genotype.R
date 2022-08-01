#' Extract SNP genotypes from aligned bamfiles using ANGSD.
#'
#' Uses the ANGSD (Analysis of Next Gen Sequencing Data) software to genotype
#' individuals across a set of bamfiles using several different filtering
#' options or genotyping methods. This function is a \emph{wrapper for ANGSD}.
#'
#' @param bamfiles character. Vector of filepaths to the bamfiles containing
#'   data for the individuals to genotype.
#' @param outfile_prefix character, default "genotypes". Prefix for the
#'   outfiles.
#' @param minInd numeric, default \code{floor(length(bamfiles)/2)}. Minimum
#'   number of individuals a locus must be sequenced in in order to call
#'   genotypes.
#' @param genotyper character, default "GATK". Name of the genotyper to use.
#'   Options: \itemize{\item{SAMtools} \item{GATK} \item{SOAPsnp} \item{phys}
#'   \item{sample: } randomly sample reads to make the genotypes (therefore
#'   posteriors are either 0, .5, or 1).} See ANGSD documentation for details.
#' @param SNP_pval numeric, default 0.00000001. Minimum p-value that must be met
#'   \emph{for the probability that a given locus is a SNP} for a locus to be
#'   genotyped.
#' @param doGeno character, default "NN". Format in which genotypes will be
#'   printed. Options: \itemize{\item{major_minor: } Print the major and minor
#'   alleles. \item{numeric: } Print alleles as 0, 1, 2, or -1 for the
#'   homozygous major, hetoerzygote, homozygous minor, or missing data,
#'   respectively. \item{NN: } Print alleles as paired nucleotides (AA, AC, CC,
#'   for example). NN denotes missing data. \item{all_posteriors: } Prints
#'   posterior probabilities for \emph{all three} genotypes at each locus.
#'   \item{called_posterior: } Prints the posterior probability for \emph{only
#'   the called genotype}. \item{binary: } Prints the posterior probabilities
#'   for the three genotypes \emph{in binary}.}
#' @param postCutoff numeric, default 0.95. The minimum posterior probability
#'   required for the most likely genotype for the genotype to be printed at for
#'   a given individual at a given locus.
#' @param minQ numeric, default 20. The minimum Phred \emph{sequencing} quality
#'   score needed for a given sequenced base on a single read to be considered
#'   during genotyping. The default, 20, corresponds to 99\% accuracy.
#' @param minMapQ numeric, default 20. The minimum Phred \emph{mapping} quality
#'   score needed for a given sequenced base on a single read to be considered
#'   during genotyping. The default, 20, corresponds to 99\% accuracy.
#' @param unzip logical, default FALSE. If TRUE, resulting genotype file will be
#'   unzipped.
#' @param doVcf logical, default FALSE. If TRUE, create a vcf file. Not
#'   supported for \code{doGeno = 'major_minor'}. If \code{doGeno} is
#'   \code{numeric} or \code{NN}, done using an 'in-house' perl script.
#'   Otherwise done using the \code{-doBcf} option in ANGSD, then converted to a
#'   .vcf using \code{bcftools}.
#' @param filter_paralogs logical, default FALSE. If true \code{ngsParalogs}
#'   will be used to detect and remove paralogs prior to genotype calling.
#' @param filter_paralogs_alpha numeric, default 0.0016. Sites with p-values
#'   below this will be considered paralogs if \code{filter_paralogs} is TRUE.
#'   The default corresponds to a likelihood ratio score of roughly 10.
#' @param filter_paralogs_buffer numeric, default 1000. Sites within this many
#'   bp of and candidate paralog sites will be excluded during genotyping if
#'   \code{filter_paralogs} is TRUE.
#' @param filter_paralogs_populations character, default rep("only_pop",
#'   length(bamfiles)). Vector of population names for each input bamfile for
#'   use during paralog filtering. Loci will be considered candidate paralogs if
#'   they are below \code{filter_paralogs_alpha} in any population.
#' @param filter_paralogs_reference character, default FALSE. If
#'   \code{filter_paralogs} is TRUE, a reference genome is needed for the
#'   creation of the target region .rf file needed by \code{SAMtools}.
#' @param rf logical, default FALSE. If TRUE, a path target region (.rf) file
#'   containing genomic regions where loci will be sequenced. An rf file can be
#'   supplied while still doing paralog filtering with if \code{filter_paralogs}
#'   is TRUE, in which case only loci from both non-paralog regions and noted in
#'   the rf file will be genotyped. Format: one region per line, noted by
#'   'chrname:start-end'. The ending position can be omitted if the rest of the
#'   remaining chromosome is a target region ('chr1:4124-' would keep every site
#'   after position 4123). A single base can be specified using 'chrname:bp-bp'
#'   (for example, 'chr1:4124-4124 would sequence only position 4124 on chr1).
#'   Multiple regions per chromosome are acceptable. An entire chromosome can be
#'   targeted by either omitting both the start and end or by noting both the
#'   first and last bp ('chr1:1-', 'chr1:', and 'chr1:1-10000' would all include
#'   the entire chr1 if it is 10000 bp long).
#' @param par numeric, default 1. Number of cores to allow ANGSD to use for
#'   genotyping.
#'   
#' @return Generates genotype files in the requested format named
#'   outfile_prefix.geno(.gz) in the directory of the first bamfile. The file
#'   will be zipped (.gz) if \code{unzip} is not TRUE. In R, returns file path
#'   to genotype output.
#'   
#' @author William Hemstrom
#' @author Michael Miller
#'
#' @export
genotype_bams <- function(bamfiles,
                          outfile_prefix = "genotypes",
                          minInd = floor(length(bamfiles)/2),
                          genotyper = "GATK",
                          SNP_pval = 0.00000001,
                          doGeno = "NN",
                          postCutoff = 0.95,
                          minQ = 20,
                          minMapQ = 20,
                          unzip = FALSE,
                          doVcf = FALSE,
                          filter_paralogs = FALSE,
                          filter_paralogs_alpha = 0.0016,
                          filter_paralogs_buffer = 1000,
                          filter_paralogs_populations = rep("only_pop", length(bamfiles)),
                          filter_paralogs_reference = FALSE,
                          rf = FALSE,
                          par = 1){
  #=============sanity checks===================
  msg <- character()

  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  if(!.check_system_install("angsd")){
    msg <- c(msg, "No angsd install located on system path.\n")
  }
  if(!.check_system_install("bcftools") & doGeno %in% c("all_posteriors", "called_posterior", "binary") & doVcf){
    msg <- c(msg, "No bcftools install located on system path. Needed for genotype posterior .vcf creation.\n")
  }

  if(filter_paralogs){
    
    if(length(filter_paralogs_populations) != length(bamfiles)){
      msg <- c(msg, "The length of 'filter_paralogs_populations' must be equal to that of the 'bamfiles'.\n")
    }
    pops <- unique(filter_paralogs_populations)
    
    if(!.check_system_install("ngsParalog")){
      msg <- c(msg, "No ngsParalog install located on system path.\n")
    }
    if(!.check_system_install("samtools")){
      msg <- c(msg, "No SAMtools install located on system path.\n")
    }
    
    if(!file.exists(filter_paralogs_reference)){
      msg <- c(msg, paste0("File not found: ", filter_paralogs_reference, "\n"))
    }
    else{
      check <- .check_is_genome(filter_paralogs_reference)
      if(is.character(check)){msg <- c(msg, check)}
      
    }
  }
  
  if(!isFALSE(rf)){
    rf <- normalizePath(rf)
    if(!file.exists(rf)){
      msg <- c(msg, paste0("Cannot locate file: ", rf, "\n"))
    }
  }
  
  bamfiles <- normalizePath(bamfiles)
  bad_bams <- which(!file.exists(bamfiles))
  if(length(bad_bams) != 0){
    msg <- c(msg, paste0("Cannot locate bamfiles: ", paste0(bad_bams, collapse = ", "), "\n"))
  }
  
  if(length(msg) > 0){
    stop(msg)
  }

  # check that the bams exist and are indexed.
  if(all(file.exists(bamfiles))){
    index_files <- paste0(bamfiles, ".bai")
    indexed <- file.exists(index_files)
    if(any(!indexed)){
      cat(sum(!indexed), "bam files not indexed. Indexing with 'samtools index'.\n")

      if(!.check_system_install("samtools")){
        stop("No samtools install located on system path.\n")
      }

      for(i in 1:sum(!indexed)){
        system(paste0("samtools index ", bamfiles[!indexed][i]))
      }
    }
  }
  else{
    msg <- c(msg, paste0("Some bamfiles not located: ",
                         paste0(bamfiles[!file.exists(bamfiles)], collapse = ", "),
                         ".\n"))
  }


  # figure out genotyper code
  genotyper <- .angsd_genotyper_key(genotyper)
  if(!is.numeric(genotyper)){
    msg <- c(msg, genotyper)
  }

  # figure out doGeno
  doGeno_table <- data.frame(number = c(1, 2, 4, 8, 16, 32),
                             word = c("major_minor", "numeric", "NN", "all_posteriors",
                                      "called_posterior", "binary"))
  bad_doGeno <- !doGeno %in% doGeno_table[,2]
  if(any(bad_doGeno)){
    msg <- c(msg, paste0("doGeno option(s) ", paste0(doGeno[bad_doGeno], collapse = ", "),
                         " not available. Is this a typo?\n"))
  }
  else if(length(doGeno) != 1){
    msg <- c(msg, "Exactly one doGeno option must be provided.\n")
  }
  else{
    doGeno <- doGeno_table[match(doGeno, doGeno_table[,2]),1]
    doGeno <- sum(doGeno)
  }

  if(doVcf){
    if(doGeno[1] == 8 | doGeno[1] == 16 | doGeno[1] == 32){
      angsd_doVcf <- TRUE
      local_doVcf <- FALSE
    }
    else if (doGeno == 1) {
      msg <- c(msg, "doVcf requires a doGeno option other than 'major_minor'.\n")
    }
    else{
      local_doVcf <- TRUE
      angsd_doVcf <- FALSE
      if(isFALSE(unzip)){
        msg <- c(msg, "'unzip' must be TRUE if a .vcf file is to be generated from a numeric or NN output.\n")
      }
      if(!.check_system_install("perl")){
        msg <- c(msg, "No perl install located on system path. This is needed if a .vcf file is to be generated from a numeric or NN output.\n")
      }
      
      doGeno <- doGeno + 1
    }
  }
  else{
    angsd_doVcf <- FALSE
    local_doVcf <- FALSE
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  

  dir <- dirname(bamfiles[1])
  outfile <- file.path(dir, outfile_prefix)

  #==============prepare to run=================
  old.scipen <- options("scipen")
  options(scipen = 999)
  script <- .fetch_a_script("angsd_genotypes.sh", "shell")
  
  # filter paralogs if requested
  if(filter_paralogs){
    rf <- .filter_paralogs(bamfiles = bamfiles, 
                           populations = filter_paralogs_populations,
                           outfile = "selected_clean_regions.txt",
                           buffer = filter_paralogs_buffer, 
                           alpha = filter_paralogs_alpha,
                           reference = filter_paralogs_reference, 
                           genotyper = genotyper, 
                           SNP_pval = SNP_pval, 
                           minQ = minQ, 
                           minMapQ = minMapQ,
                           cleanup = TRUE,
                           rf = rf,
                           par = par)
  }

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
                  ifelse(angsd_doVcf, 1, 0),
                  ifelse(file.exists(rf), 1, 0),
                  ifelse(file.exists(rf), rf, "NA"),
                  par), collapse = " ")
                )

  system(cmd)
  
  # if returning a vcf, convert and do so
  if(angsd_doVcf){
    system(paste0("bcftools convert -O v -o ", outfile, ".vcf ", outfile, ".bcf"))
    output_filename <- paste0(outfile, ".vcf")
    return(output_filename)
  }
  
  output_filename <- paste0(outfile, ".geno.gz")

  if(unzip){
    system(paste0("gunzip ", outfile, ".geno.gz"))
    output_filename <- paste0(outfile, ".geno")
  }
  
  if(local_doVcf){ # only TRUE if unzip is TRUE due to sanity checks, so don't need to worry about a zipped input.
    vcf_script <- .fetch_a_script("ConvertGenosToVCF.pl", "perl")
    tf <- tempfile()
    write.table(basename(bamfiles), tf, col.names = FALSE, row.names = FALSE, quote = FALSE)
    system(paste0("perl ", vcf_script, " ",
                  paste0(outfile, ".geno"), " ",
                  paste0(outfile, ".vcf"), " ",
                  ifelse(doGeno == 5, "NN ", "numeric "),
                  tf))
    file.remove(tf)
    return(paste0(outfile, ".vcf"))
  }
  
  options(scipen = old.scipen)
  
  return(output_filename)
}
