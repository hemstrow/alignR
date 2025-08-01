#' @param read_metadata a data.frame or equivalent with a row for each \code{RA}
#'   read file and a columns, in order, indicating \itemize{
#'   \item{library}
#'   \item{sampleID/name}
#'   \item{sequencing platform}
#'   \item{flowcell}
#'   \item{sequencing lane}}
#'   Column names are optional but recommended for organizational purposes.
#'   The 'flowcell' and 'sequencing lane' columns can be left as NA,
#'   in which case they will be guessed from the fastq header. This may create
#'   errors if fastq headers are non-standard (flow cell and lane are expected
#'   in the third and fourth slots, respectively (separated by ':')).
run_rmake_pipeline <- function(RA, RB = NULL,
                               read_metadata,
                               reference,
                               chr_size_cutoff,
                               max_chunk_size,
                               merge_list = NULL,
                               until = NULL,
                               downsample_args = list(low_cut = 0, high_cut = Inf),
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
                                                 platform = "ILLUMINA"),
                               HaplotypeCaller_args = list(mem = 2),
                               GenomicsDBImport_args = list(batch_size = 5),
                               VariantFiltration_args = list(QD = 2,
                                                             FS = 60,
                                                             SOR = 3,
                                                             MQ = 40,
                                                             MQRankSum = -12.5,
                                                             ReadPosRankSum = -8,
                                                             min_genotype_quality = 13),
                               slurm_profile = list(genome_index = NULL,
                                                    trim = NULL,
                                                    align = NULL,
                                                    merge = NULL,
                                                    downsample = NULL)){

  # rmake needs to know where the script is that generates the makefile; so save a dummy copy of this funciton first
  sink("Makefile.R")
  print.function(run_rmake_pipeline)
  sink()


  #=======================sanity checks===========================================
  slurm_profile_needs <- c("genome_index", "trim", "align", "merge", "downsample", "HaplotypeCaller")
  for(i in 1:length(slurm_profile_needs)){
    if(!slurm_profile_needs[i] %in% names(slurm_profile)){
      nl <- list(temp = NULL)
      names(nl) <- slurm_profile_needs[i]
      slurm_profile <- c(slurm_profile, nl)
    }
  }
  slurm_profile <- slurm_profile[which(names(slurm_profile) %in% slurm_profile_needs)]

  sample_names <- rmake::replaceSuffix(RA, "")
  sample_namesB <- rmake::replaceSuffix(RB, "")



  #=======================generate scatter info===================================
  # Computationally light, this is the only real "computation" step done prior to the actual snakemake. Specifically, this will
  # 1) faidx the genome only if it hasn't been already;
  # 2) Generate scatter metadata from the .fai data for use in the rmake to follow; and
  # 3) write scatter lists ONLY if they don't already exist or are different (implying a change in chunking or ref genome)


  # first write region GATK-style .list files for HaplotypeCaller (https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists)
  scatter_info <- prep_scatter_intervals(reference, chr_size_cutoff, max_chunk_size)
  scatter_info <- data.table::as.data.table(scatter_info)

  # figure out what list files we need and re-write ONLY if they don't already exist and are identical to what we would make here. This allows this function to be re-run and work like a proper make run.
  scatter_groups <- unique(scatter_info$scatter_idx)
  for(i in 1:length(scatter_groups)){
    need_write <- FALSE
    tfile <- file.path(dirname(RA[1]), paste0(scatter_groups[i], ".list"))
    twrite <- unlist(scatter_info[scatter_idx == scatter_groups[i], list(paste0(chrom, ":", start, "-", end))])


    if(file.exists(tfile)){
      comp <- readLines(tfile)
      if(!all(comp == twrite)){
        need_write <- TRUE
      }
    }
    else{
      need_write <- TRUE
    }

    if(need_write){writeLines(twrite, tfile)}
  }


  #=======================pre-processing rules===================================================
  # key notes: sample names will hold sample names (B holds R2 reads while that is used)
  #                These may be updated following merging
  #            RA (and RB at the start) will hold current output file names, including extensions.

  # prep genome
  ref_bn <- gsub("\\.gz$", "", reference)
  ref_bn <- rmake::replaceSuffix(ref_bn, "")
  gprep_script <- .fetch_a_script("gprep.R", "rmake")
  gprep_rule <- rmake::rRule(target = c(paste0(reference, c(".sa",
                                                          ".amb",
                                                          ".bwt",
                                                          ".ann",
                                                          ".bwt",
                                                          ".fai",
                                                          ".pac")),
                                        paste0(ref_bn, ".dict")),
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

  trim_rule <- rmake::expandTemplate(trim_rule, data.frame(RA = sample_names,
                                                           RB = sample_namesB))
  RA <- rmake::replaceSuffix(RA, "_trimmed.fastq")
  RB <- rmake::replaceSuffix(RB, "_trimmed.fastq")

  # align
  browser()
  align_script <- .fetch_a_script("align.R", "rmake")

  align_rule <- rmake::rRule(target = c("$[out].sort.flt.bam",
                                        "$[out].sort.flt.bam.bai",
                                        "$[out].sort.flt.bam.flagstat",
                                        "$[out].sort.flt.bam.stats",
                                        "$[out].markdup.bam.flagstat"),
                        script = align_script,
                        depends = c("$[RA]", "$[RB]", paste0(reference, c(".sa",
                                                                          ".amb",
                                                                          ".bwt",
                                                                          ".ann",
                                                                          ".bwt",
                                                                          ".pac"))),
                        params = .fill_args_defaults(c(list(reference = reference,
                                                            read_metadata = cbind(sn = RA, read_metadata)),
                                                       c(align_args,
                                                         slurm_profile$align)),
                                                     func = align_reference,
                                                     skip = c("RA_fastqs", "RB_fastqs")))
  align_rule <- rmake::expandTemplate(align_rule, data.frame(out = sample_names,
                                                             RA = RA,
                                                             RB = RB))
  RA <- paste0(sample_names, ".sort.flt.bam")

  # merge anything that needs merging
  if(is.list(merge_list)){
    merge_script <- .fetch_a_script("merge.R", "rmake")
    merge_rule <- vector("list", length(merge_list))



    # can't use expandTemplate to make these rules, since it takes multiple outputs but only one input
    for(i in 1:length(merge_rule)){

      # fix merge list with up-to-date prefixes
      tsamps <- sample_names[which(sample_names %in% rmake::replaceSuffix(merge_list[[i]], ""))]
      tml <- list(paste0(tsamps, ".sort.flt.bam"))
      names(tml) <- names(merge_list)[i]

      # make rule

      merge_rule[[i]] <- rmake::rRule(target = c(file.path(dirname(RA)[1], paste0(names(merge_list)[i], ".bam")),
                                                 file.path(dirname(RA)[1], paste0(names(merge_list)[i], ".bam.bai"))),
                                      script = merge_script,
                                      depends = c(paste0(tsamps, ".sort.flt.bam"),
                                                  paste0(tsamps, ".sort.flt.bam.bai")),
                                      params = .fill_args_defaults(args = c(list(file_list = tml),
                                                                            slurm_profile$merge),
                                                                   func = merge_bams))
    }

    # update names
    unmerged <- which(!sample_names %in% rmake::replaceSuffix(unlist(merge_list), ""))
    RA <- c(file.path(dirname(RA)[1], paste0(names(merge_list), ".bam")), RA[unmerged])
    sample_names <- c(file.path(dirname(RA)[1], names(merge_list)), sample_names[unmerged])

  }
  else{
    merge_rule <- NULL
  }

  # downsample
  if(downsample_args$low_cut > 0 | downsample_args$high_cut < Inf){
    downsample_script <- .fetch_a_script("downsample.R", "rmake")
    downsample_rule <- rmake::rRule(target = c("$[samples].sub.bam",
                                               "$[samples].sub.bam.bai"),
                                    script = downsample_script,
                                    depends = c("$[RA]", "$[RA].bai"),
                                    params = .fill_args_defaults(c(list(reference = reference),
                                                                   c(downsample_args,
                                                                     slurm_profile$downsample)),
                                                                 func = downsample_bams,
                                                                 skip = "bams"))
    downsample_rule <- rmake::expandTemplate(downsample_rule, data.frame(RA = RA,
                                                                         samples = sample_names))

    RA <- paste0(sample_names, ".sub.bam")
  }
  else{
    downsample_rule <- NULL
  }

  # the next few steps are scattered
  # basic idea is: scatter chromosomes and scaffold groups, run Haplotype caller for each interval, import into shared dbs for scaffold groups and chromosomes, then genotype for each interval

  #=======================GATK multistep==========================================
  # Haplotype caller
  browser()
  HC_script <- .fetch_a_script("HaplotypeCaller.R", "rmake")

  HC_rule <- rmake::rRule(target = c("$[sg]-$[rf].hapcalls.gvcf.gz",
                                     "$[sg]-$[rf].hapcalls.gvcf.gz.tbi"),
                          script = HC_script,
                          depends = c("$[RA]", file.path(dirname(RA[1]), "$[rf].list"), "$[RA].bai",
                                      paste0(reference, ".fai"),
                                      paste0(ref_bn, ".dict")),
                          params = .fill_args_defaults(c(list(reference = reference),
                                                         slurm_profile$HaplotypeCaller,
                                                         HaplotypeCaller_args),
                                                       func = run_HaplotypeCaller,
                                                       skip = c("bamfiles", "region")))

  HC_fill_df <- merge(data.frame(sg = sample_names, RA = RA),
                            data.frame(rf = scatter_groups))
  HC_rule <- rmake::expandTemplate(HC_rule, vars = HC_fill_df)

  rmake::makefile(c(list(gprep_rule), trim_rule, align_rule, merge_rule, downsample_rule, HC_rule), fileName = "Makefile", makeScript = "Makefile.R")
  rmake::make()



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
