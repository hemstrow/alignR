run_rmake_pipeline <- function(RA, RB = NULL,
                               reference,
                               merge_list = NULL,
                               min_chr_size = 0,
                               chunk_size = "chr",
                               until = NULL,
                               min_coverage = 0,
                               downsample = FALSE,
                               fastp_args = list(adapter_r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                                                 adapter_r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                                                 cut_front = FALSE,
                                                 cut_front_mean_quality = NULL,
                                                 cut_tail = FALSE,
                                                 cut_tail_mean_quality = NULL,
                                                 cut_right = 4,
                                                 cut_right_mean_quality = 20,
                                                 trim_poly_g = 10,
                                                 trim_poly_x = FALSE),
                               align_args = list(mapQ = 5,
                                                 remove_duplicates = TRUE,
                                                 remove_improper_pairs = TRUE,
                                                 par = 1),
                               add_RGS_args = list(platform = "ILLUMINA"),
                               GenomicsDBImport_args = list(batch_size = 5),
                               VariantFiltration_args = list(QD = 2,
                                                             FS = 60,
                                                             SOR = 3,
                                                             MQ = 40,
                                                             MQRankSum = -12.5,
                                                             ReadPosRankSum = -8,
                                                             min_genotype_quality = 13),
                               slurm_profile = NULL){

  # rmake needs to know where the script is that generates the makefile; so save a dummy copy of this funciton first
  sink("Makefile.R")
  print.function(run_rmake_pipeline)
  sink()

  #=======================rules===================================================
  # prep genome
  gprep_script <- .fetch_a_script("gprep.R", "rmake")
  gprep_rule <- rmake::rRule(target = paste0(reference, c(".sa",
                                                          ".amb",
                                                          ".bwt",
                                                          ".ann",
                                                          ".bwt",
                                                          ".fai",
                                                          ".pac",
                                                          ".len")),
                             script = gprep_script,
                             depends = reference,
                             params = c(list(), slurm_profile$genome_index))
  # trimming
  trim_script <- .fetch_a_script("trim.R", "rmake")

  trim_rule <- rmake::rRule(target = c("$[RA]_trimmed.fastq", "$[RB]_trimmed.fastq"),
                            script = trim_script,
                            depends = c("$[RA].fastq", "$[RB].fastq"),
                            params = .fill_args_defaults(c(c(fastp_args,
                                                             slurm_profile$trim)),
                                                         func = trim_fastp,
                                                         skip = c("RA_fastqs", "RB_fastqs")))

  trim_rule <- rmake::expandTemplate(trim_rule, data.frame(RA = rmake::replaceSuffix(RA, ""),
                                                           RB = rmake::replaceSuffix(RB, "")))

  # alignment
  align_script <- .fetch_a_script("align.R", "rmake")

  align_rule <- rmake::rRule(target = c("$[RA].sort.flt.bam",
                                        "$[RA].sort.flt.bam.bai",
                                        "$[RA].sort.flt.bam.flagstat",
                                        "$[RA].sort.flt.bam.stats",
                                        "$[RA].markdup.bam.flagstat"),
                        script = align_script,
                        depends = c("$[RA]_trimmed.fastq", "$[RB]_trimmed.fastq", rmake::targets(gprep_rule)),
                        params = .fill_args_defaults(c(list(reference = reference),
                                                       c(align_args,
                                                         slurm_profile$align)),
                                                     func = align_reference,
                                                     skip = c("RA_fastqs", "RB_fastqs")))
  align_rule <- rmake::expandTemplate(align_rule, data.frame(RA = rmake::replaceSuffix(RA, ""),
                                                             RB = rmake::replaceSuffix(RB, "")))

  rmake::makefile(c(list(gprep_rule), trim_rule, align_rule), fileName = "Makefile", makeScript = "Makefile.R")
  rmake::make()

  # downsample
  downsample_script <- .fetch_a_script("downsample.R", "rmake")
  downsmple_rule <- rmake::rRule(target = c("$[RA].sort.flt.sub.bam",
                                            "$[RA].sort.flt.sub.bam.bai"),
                                 script = downsample_script,
                                 depends = c("$[RA]_trimmed.fastq", "$[RB]_trimmed.fastq", rmake::targets(gprep_rule)),
                                 params = .fill_args_defaults(c(list(reference = reference),
                                                                c(align_args,
                                                                  slurm_profile$align)),
                                                              func = align_reference,
                                                              skip = c("RA_fastqs", "RB_fastqs")))

  # scattering

  # GATK--multistep.

  #=======================makefile function=======================================
  job <- list(trim, align)

  makefile(job, "Makefile")


}

# fill in default args from a list of args for a function, optionally appending
# any args starting in "slurm_"
.fill_args_defaults <- function(args, func, append_slurm = TRUE, skip = NULL){
  defaults <- formals(func)
  for(i in 1:length(args)){
    if(names(args)[i] %in% names(defaults)){
      defaults[[names(args)[i]]] <- args[[i]]
    }
  }

  # combine them to add slurm args, checking they are legit first
  if(append_slurm){
    missing <- which(names(args) %in% names(defaults))
    missing <- args[missing]
    if(length(missing) > 0){
      OKmissing <- missing[grep("^slurm_", missing)]
      defaults <- c(defaults, OKmissing)
    }
  }

  # if skipping any, these args are omitted (they are processed some other way)
  if(!is.null(skip)){
    if(any(names(defaults) %in% skip)){
      defaults <- defaults[-which(names(defaults) %in% skip)]
    }
  }
  return(defaults)
}

# calling function that runs a function provided with args, some of which might be slurm args.
# this should run within the rule scripts sourced by rmake::rRule()
.rule_caller_w_slurm <- function(func, params, jobname = NULL){
  # grab sets of params
  slurm_params <- params[grep("^slurm_", params)]
  if(length(slurm_params) > 0){
    params <- params[-grep("^slurm_", params)]
  }

  if(length(slurm_params) == 0){
    params <- params[which(names(params) %in% names(formals(func)))]
    do.call(func, args = params)
  }

  # call with rslurm
  else{
    params <- params[which(names(params) %in% names(formals(func)))]
    rslurm::slurm_call(func, params = params, jobname = jobname, slurm_options = slurm_params)
  }
}

.console_hline <- function(char = "="){
  return(paste0(rep(char, getOption("width")), collapse = ""))
}
