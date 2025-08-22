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
                               HaplotypeCaller_args = list(mem = 2,
                                                           min_base_quality_score = 33,
                                                           minimum_mapping_quality = 20),
                               GenomicsDBImport_args = list(mem = 4,
                                                            batch_size = 5),
                               GenotypeGVCFs_args = list(mem = 4,
                                                         max_alternate_alleles = 2,
                                                         new_qual = TRUE),
                               VariantFiltration_args = list(mem = 2,
                                                             QD = 2,
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

  # rmake needs to know where the script is that generates the makefile; so save a dummy copy of this function alongside called args.
  env <- c(as.list(environment()))
  sink("Makefile.R")
  cat("# alignR run_rmake_pipeline call.\n#=====================================\n#Arguments:\n")
  for(i in 1:length(env)){
    cat(paste0("# ", names(env)[i], "\n#\t\t", paste0(unlist(env[[i]]), collapse = ", "), "\n"))
  }
  cat("# alignR run_rmake_pipeline call.\n#=====================================\n#Function:\n")
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

  if(!.check_system_install("gatk4")){
    stop("No gatk4 install located on system path.\n")
  }
  if(!.check_system_install("picard")){
    stop("No picard install located on system path.\n")
  }
  if(!.check_system_install("samtools")){
    stop("No SAMtools install located on system path.\n")
  }
  if(!.check_system_install("bcftools")){
    stop("No bcftools install located on system path.\n")
  }
  if(!.check_system_install("bwa")){
    stop("No bwa install located on system path.\n")
  }
  if(!.check_system_install("vcftools")){
    stop("No vcftools install located on system path.\n")
  }

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
  # HaplotypeCaller
  HC_script <- .fetch_a_script("HaplotypeCaller.R", "rmake")

  HC_rule <- rmake::rRule(target = c("$[sg]-$[rf].hapcalls.gvcf.gz",
                                     "$[sg]-$[rf].hapcalls.gvcf.gz.tbi"),
                          script = HC_script,
                          depends = c("$[RA]", file.path(dirname(RA[1]), "$[rf].list"), "$[RA].bai",
                                      reference,
                                      paste0(ref_bn, ".dict")),
                          params = .fill_args_defaults(c(list(reference = reference),
                                                         slurm_profile$HaplotypeCaller,
                                                         HaplotypeCaller_args),
                                                       func = run_HaplotypeCaller,
                                                       skip = c("bamfiles", "region")))

  HC_fill_df <- merge(data.frame(sg = sample_names, RA = RA),
                            data.frame(rf = scatter_groups))
  HC_rule <- rmake::expandTemplate(HC_rule, vars = HC_fill_df)

  # CenomicsDBImport
  # this runs on each scaffold group/chromosome.
  GDBI_script <- .fetch_a_script("GenomicsDBImport.R", "rmake")
  ## rule is a bit complicated, so made with for loop(s). Note also that this doesn't use the genomicsDBImport function in alignR because that uses hapmap files that we don't have.
  GDBI_rule <- vector("list", 0)
  gdbis <- character()
  if(any(scatter_info$type == "chrom")){
    uchroms <- unique(scatter_info[type == "chrom",]$chrom)
    for(i in 1:length(uchroms)){
      tscats <- scatter_info[type == "chrom" & chrom == uchroms[i],]$scatter_idx
      tdeps <- expand.grid(sample_names, tscats)
      tdeps <- paste0(tdeps[[1]], "-", tdeps[[2]], ".hapcalls.gvcf.gz")
      GDBI_rule[i] <- list(rmake::rRule(target = file.path(dirname(RA[1]), paste0(uchroms[i], "_db")),
                                        script = GDBI_script,
                                        depends = c(tdeps, paste0(tdeps, ".tbi"),
                                                    reference,
                                                    paste0(ref_bn, ".dict")),
                                        params = c(list(reference = reference, L = uchroms[i]),
                                                   slurm_profile$GenomicsDBImport,
                                                   GenomicsDBImport_args)))
      gdbis <- c(gdbis, rmake::targets(GDBI_rule[[i]]))
    }
  }
  if(any(scatter_info$type == "scaf")){
    sgs <- unique(scatter_info[type == "scaf",]$scatter_idx)
    for(i in 1:length(sgs)){
      tdeps <- expand.grid(sample_names, sgs[i])
      tdeps <- paste0(tdeps[[1]], "-", tdeps[[2]], ".hapcalls.gvcf.gz")
      GDBI_rule <- c(GDBI_rule, list(rmake::rRule(target = file.path(dirname(RA[1]), paste0(sgs[i], "_db")),
                                                  script = GDBI_script,
                                                  depends = c(tdeps, paste0(tdeps, ".tbi"),
                                                              reference,
                                                              paste0(ref_bn, ".dict")),
                                                  params = c(list(reference = reference,
                                                                  L = file.path(dirname(RA[1]),
                                                                                paste0(scatter_groups, ".list"))),
                                                             slurm_profile$GenomicsDBImport,
                                                             GenomicsDBImport_args))))

      gdbis <- c(gdbis, rmake::targets(GDBI_rule[[length(GDBI_rule)]]))
    }
  }

  # GenotypeGVCFs
  # as before this uses a slightly different syntax (no hapmaps) than the alignR wrapper function, so this gets it's own version.
  genotype_script <- .fetch_a_script("GenotypeGVCFs.R", "rmake")
  genotype_rule <- list()
  if(any(scatter_info$type == "chrom")){
    chrom_geno_rule <- rmake::rRule(target = file.path(dirname(RA[1]), "$[scatter_idx].vcf.gz"),
                                    script = genotype_script,
                                    depends = c(file.path(dirname(RA[1]), "$[scatter_idx].list"),
                                                file.path(dirname(RA[1]), "$[chrom]_db"),
                                                reference,
                                                paste0(ref_bn, ".dict")),
                                    params = c(slurm_profile$GenotypeGVCFs,
                                               GenotypeGVCFs_args))

    genotype_rule <- c(genotype_rule, rmake::expandTemplate(chrom_geno_rule, vars = as.data.frame(scatter_info[type == "chrom", .(scatter_idx, chrom)])))
  }
  if(any(scatter_info$type == "scaf")){
    scat_geno_rule <- rmake::rRule(target = file.path(dirname(RA[1]), "$[scatter_idx].vcf.gz"),
                                    script = genotype_script,
                                    depends = c(file.path(dirname(RA[1]), "$[scatter_idx].list"),
                                                file.path(dirname(RA[1]), "$[scatter_idx]_db"),
                                                reference,
                                                paste0(ref_bn, ".dict")),
                                    params = c(slurm_profile$GenotypeGVCFs,
                                               GenotypeGVCFs_args))

    genotype_rule <- c(genotype_rule, rmake::expandTemplate(scat_geno_rule, vars = as.data.frame(scatter_info[type == "scaff",]$scatter_idx)))
  }

  # Variant Filtration
  ## we can use the pre-existing function here, no problems
  browser()
  VF_script <- .fetch_a_script("VariantFiltration.R", "rmake")

  VF_rule <- rmake::rRule(target = c("$[scat]_hard_filt_pass.recode.vcf"),
                          script = VF_script,
                          depends = c("$[scat].vcf.gz",
                                      reference,
                                      paste0(ref_bn, ".dict")),
                          params = .fill_args_defaults(c(slurm_profile$VariantFiltration,
                                                         VariantFiltration_args),
                                                       func = run_VariantFiltration,
                                                       skip = c("vcfs", "reference")))
  VF_rule <- rmake::expandTemplate(VF_rule, vars = data.frame(scat = file.path(dirname(RA[1]), scatter_groups)))



  rmake::makefile(c(list(gprep_rule), trim_rule, align_rule, merge_rule, downsample_rule,
                    HC_rule, GDBI_rule, genotype_rule, VF_rule),
                  fileName = "Makefile", makeScript = "Makefile.R")
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
