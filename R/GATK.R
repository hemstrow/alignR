#' @describeIn genotype_bams_GATK Build ".gvcf" haplotype files for each sample
#'   with \code{GATK}'s \code{HaplotypeCaller}.
#' @export
run_HaplotypeCaller <- function(bamfiles, reference, mem, par = 1, java_path = "java", gatk4_path = "gatk.jar"){
  #==========sanity checks=========
  msg <- character()
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

  if(length(msg) > 0){
    stop(msg)
  }

  #=======================run=========
  script <- .fetch_a_script("run_HaplotypeCaller.sh", "shell")


  temp_dir <- tempdir()


  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # divide up into ncore chunks
  iters <- length(bamfiles)
  it_par <- (1:iters)%%par
  chunks <- split(bamfiles, it_par)

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    valid <- data.frame(bam = chunks[[q]], valid = character(length(chunks[[q]])), msg = character(length(chunks[[q]])))
    for(i in 1:length(chunks[[q]])){
      # run Haplotype Caller

      ttempdir <- file.path(temp_dir, paste0("HC_q", q, "_i", i))
      dir.create(ttempdir)
      cmd <- paste0("bash ", script, " ",
                    chunks[[q]][i], " ",
                    reference, " ",
                    ttempdir, " ",
                    mem, " ",
                    java_path, " ",
                    gatk4_path, " ")

      system(cmd)
      unlink(ttempdir, recursive = TRUE)

      # validate
      cmd <- paste0(java_path, " -jar ", gatk4_path, " ValidateVariants ",
                   "-R ", reference, " ",
                   "-gvcf TRUE ",
                   "--verbosity ERROR ",
                   "-V ", paste0(chunks[[q]][i], ".hapcalls.gvcf.gz"),
                   " 2>&1")

      valid[i,3] <- paste0(system(cmd, intern = TRUE), collapse = "\n")
      valid$valid[i] <- !grepl("ERR", valid[i,3])
    }

    valid
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  # check validity
  valid <- data.table::rbindlist(output)
  valid[,2] <- as.logical(valid[,2])

  if(any(!valid[,2])){
    warning(paste0("Some HaplotypeCaller runs produced invalid vcf outputs. This can occur for a few reasons, but too small memory limits are a common cause. Invalid outputs:",
                   paste0(valid[!valid[,2],1], sep = "\n\t"),
                   "\nMessages for bad runs dumped to: ./ValidateVariants.err"))
    writeLines(valid[!valid[,2], 3], "./ValidateVariants.err", sep = "\n\n===========================================================\n\n")
  }


  # save a hapmap
  sample_map <- data.frame(samp = basename(tools::file_path_sans_ext(bamfiles)), file = normalizePath(paste0(bamfiles, ".hapcalls.gvcf.gz")))
  data.table::fwrite(sample_map, paste0(bamfiles[1], "_hapmap.txt"), sep = "\t", col.names = F, row.names = F)
  cat("hapmap saved to: ", paste0(bamfiles[1], "_hapmap.txt"), "\n")

  # return info
  return(paste0(bamfiles[1], "_hapmap.txt"))
}

#' @describeIn genotype_bams_GATK Add read groups to bam files with
#'   \code{picard}.
#' @export
add_RGS <- function(bamfiles, fastqs, par = 1, platform = "ILLUMINA", java_path = "java", picard_path = "picard.jar"){
  #==========sanity checks=========
  msg <- character()
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

  # check that the fastqs all exist
  fastqs <- normalizePath(fastqs)
  bad_fastqs <- which(!file.exists(fastqs))
  if(length(bad_fastqs) != 0){
    msg <- c(msg, paste0("Cannot locate fastqs: ", paste0(bad_fastqs, collapse = ", "), "\n"))
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #=====================run============

  script <- .fetch_a_script("add_RGs.sh", "shell")

  iters <- length(bamfiles)
  it_par <- (1:iters)%%par
  chunks <- split(bamfiles, it_par)
  chunks_f <- split(fastqs, it_par)

  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      cmd <- paste0("bash ", script, " ",
                    chunks[[q]][i], " ",
                    chunks_f[[q]][i], " ",
                    basename(tools::file_path_sans_ext(chunks[[q]][i])), " ",
                    java_path, " ",
                    picard_path, " ",
                    platform)

      system(cmd)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(paste0(tools::file_path_sans_ext(bamfiles), ".RG.bam"))
}

#' Make .bed files for chromosomes or genomic chunks.
#'
#' Generate ".bed" files describing eithe chromosomes or for automatically
#' tabulated chunks of the genome of a given size.
#'
#' @param reference Character. Path to reference genome.
#' @param min_chr_size Numeric, default 0. Minimum size for chromosomes/scaffolds.
#'   Scaffolds smaller than this will be skipped when writing ".bed" files or
#'   chunking the genome.
#' @param chunk_size Numeric or "chr", default "chr". Size of chunks by which
#'   to split and write ".bed" files, in base pairs. If "chr", a ".bed" file
#'   will be written for each chromosome/scaffold.
#' @param outdir Character, default ".". Path to directory where ".bed" files
#'   will be written.
#'
#' @export
#' @author William Hemstrom
#' @return A vector of written ".bed" file paths.
make_region_beds <- function(reference, min_chr_size = 0, chunk_size = "chr", outdir = "."){

  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  # figure out chr info:
  if(grepl("\\.gz$", reference)){
    cmd <- paste0("zcat ", reference)
  }
  else{
    cmd <- paste0("cat ", reference)
  }

  cmd <- paste0(cmd,
                r"( | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chr_lengths.txt)")
  system(cmd)

  d <- readLines("chr_lengths.txt")
  if(any(d == "")){
    d <- d[-which(d == "")]
  }

  dl <- strsplit(d, "\t")
  dl <- unlist(purrr::map(dl, 2))


  bed <- data.frame(gsub(" .+", "", gsub("\t.+", "", d)), 0, as.numeric(dl))

  if(any(bed[,3] < min_chr_size)){
    bed <- bed[-which(bed[,3] < min_chr_size),]
  }

  .batch_genotyping <- function(chr_info, chunk_size){
    old_scipen <- options("scipen")
    options(scipen = 999)

    rf_prog <- 1
    for(i in 1:nrow(chr_info)){

      # if writing a whole chr
      if(chunk_size == "chr"){
        writeLines(paste0(chr_info[i,], collapse = "\t"), file.path(outdir, paste0("region", rf_prog, ".bed")))
        rf_prog <- rf_prog + 1
        next
      }

      # otherwise chunk
      n_chuncks <- ceiling(chr_info[,3][i]/chunk_size)
      starts <- seq(0, chr_info[,3][i] - 1, length.out = n_chuncks + 1)
      starts[1:(length(starts) - 1)] <- floor(starts[1:(length(starts) - 1)])


      for(j in 1:(length(starts) - 1)){
        cmd <- paste0(chr_info[i,1], "\t", starts[j], "\t", starts[j + 1])
        writeLines(cmd, file.path(outdir, paste0("region", rf_prog, ".bed")), sep = "\n")
        rf_prog <- rf_prog + 1
      }
    }

    options(scipen = old_scipen$scipen)
    return(file.path(outdir, paste0("region", 1:(rf_prog - 1), ".bed")))
  }

  out <- .batch_genotyping(bed, chunk_size)
  return(out)
}

#' @describeIn genotype_bams_GATK Build genome databases using \code{GATK}'s
#'   \code{GenomicsDBImport}.
#' @export
run_GenomicsDBImport <- function(hapmap, bedfiles, mem, par = 1, batch_size = 5, java_path = "java", gatk4_path = "gatk.jar"){
  #============sanity checks==================
  msg <- character()

  bedfiles <- normalizePath(bedfiles)
  bad_beds <- which(!file.exists(bedfiles))
  if(length(bad_beds) != 0){
    msg <- c(msg, paste0("Cannot locate bedfiles: ", paste0(bad_beds, collapse = ", "), "\n"))
  }

  hapmap <- normalizePath(hapmap)
  if(!file.exists(hapmap)){
    msg <- c(msg, paste0("Cannot locate hapmap at path ", hapmap))
  }
  else{
    gvcfs <- fread(hapmap)[,2]
    bad_gvcfs <- which(!file.exists(gvcfs))
    if(length(bad_gvcfs) != 0){
      msg <- c(msg, paste0("Cannot locate gvcfs: ", paste0(bad_gvcfs, collapse = ", "), "\n"))
    }

    tbis <- paste0(gvcfs, ".tbi")
    bad_tbis <- which(!file.exists(tbis))
    if(length(bad_tbis) != 0){
      msg <- c(msg, paste0("Some 'tbi' index files not found for gvcfs, did these hapcalls complete?: ", paste0(bad_tbis, collapse = ", "), "\n"))
    }
  }

  if(length(msg) > 0){
    stop(msg)
  }

  #==========run==========

  script <- .fetch_a_script("run_GenomicsDBImport.sh", "shell")


  temp_dir <- tempdir()


  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # divide up into ncore chunks
  iters <- length(bedfiles)
  it_par <- (1:iters)%%par
  chunks <- split(bedfiles, it_par)

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      ttempdir <- file.path(temp_dir, paste0("GDBImport_q", q, "_i", i))
      dir.create(ttempdir)
      cmd <- paste0("bash ", script, " ",
                    hapmap, " ",
                    chunks[[q]][i], " ",
                    mem, " ",
                    ttempdir, " ",
                    java_path, " ",
                    gatk4_path, " ",
                    batch_size)

      system(cmd)
      unlink(ttempdir, recursive = TRUE)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  # prepare file list, return
  return(paste0(tools::file_path_sans_ext(bedfiles), "_db"))
}

#' @describeIn genotype_bams_GATK Generate called genotypes from GenomeDBs using
#'  \code{GATK}'s \code{GenotypeGVCFs}.
#' @export
run_GenotypeGVCFs <- function(bedfiles, reference, mem, par = 1, java_path = "java", gatk4_path = "gatk.jar"){
  #============sanity checks==================
  msg <- character()

  bedfiles <- normalizePath(bedfiles)
  bad_beds <- which(!file.exists(bedfiles))
  if(length(bad_beds) != 0){
    msg <- c(msg, paste0("Cannot locate bedfiles: ", paste0(bad_beds, collapse = ", "), "\n"))
  }

  # check that dbs exist
  genome_dbs <- paste0(tools::file_path_sans_ext(bedfiles), "_db")
  bad_dbs <- which(!dir.exists(genome_dbs))
  if(length(bad_dbs) != 0){
    msg <- c(msg, paste0("Cannot locate GenomicsDBs: ", paste0(bad_dbs, collapse = ", "), "\nDBs are expected to be named matching .bed files, ending in _db"))
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #==========run============
  script <- .fetch_a_script("run_GenotypeGVCFs.sh", "shell")

  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # divide up into ncore chunks
  iters <- length(bedfiles)
  it_par <- (1:iters)%%par
  chunks <- split(bedfiles, it_par)

  temp_dir <- tempdir()

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      # run GenotypeGVCFs

      ttempdir <- file.path(temp_dir, paste0("GenotypeGVCFs_q", q, "_i", i))
      dir.create(ttempdir)
      cmd <- paste0("bash ", script, " ",
                    chunks[[q]][i], " ",
                    reference, " ",
                    mem, " ",
                    ttempdir, " ",
                    java_path, " ",
                    gatk4_path)

      system(cmd)
      unlink(ttempdir, recursive = TRUE)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(paste0(tools::file_path_sans_ext(bedfiles), "_raw.vcf"))
}

#' @describeIn genotype_bams_GATK Filter called genotypes using \code{GATK}'s
#'   \code{VariantFiltration}.
run_VariantFiltration <- function(vcfs, reference, mem,
                                  par = 1,
                                  QD = 2,
                                  FS = 60,
                                  SOR = 3,
                                  MQ = 40,
                                  MQRankSum = -12.5,
                                  ReadPosRankSum = -8,
                                  min_genotype_quality = 13,
                                  java_path = "java", gatk4_path = "gatk.jar"){
  #============sanity checks==================
  msg <- character()

  vcfs <- normalizePath(vcfs)
  bad_vcfs <- which(!file.exists(vcfs))
  if(length(bad_vcfs) != 0){
    msg <- c(msg, paste0("Cannot locate vcfs: ", paste0(bad_vcfs, collapse = ", "), "\n"))
  }

  if(length(msg) > 0){
    stop(msg)
  }

  #==========run==========
  script <- .fetch_a_script("run_VariantFiltration.sh", "shell")

  # register parallel architecture
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)

  # divide up into ncore chunks
  iters <- length(vcfs)
  it_par <- (1:iters)%%par
  chunks <- split(vcfs, it_par)

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      # run Haplotype Caller
      cmd <- paste0("bash ", script, " ",
                    chunks[[q]][i], " ",
                    reference, " ",
                    format(QD, nsmall = 1), " ",
                    format(FS, nsmall = 1), " ",
                    format(SOR, nsmall = 1), " ",
                    format(MQ, nsmall = 1), " ",
                    format(MQRankSum, nsmall = 1), " ",
                    format(ReadPosRankSum, nsmall = 1), " ",
                    min_genotype_quality, " ",
                    mem, " ",
                    java_path, " ",
                    gatk4_path)

      system(cmd)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  return(paste0(tools::file_path_sans_ext(vcfs), "_hard_filt_pass.recode.vcf"))
}

#' @describeIn genotype_bams_GATK Prepare a reference genome for genotyping with
#'   \code{picard} and other indexers.
#' @export
prep_genome_GATK <- function(reference,
                             java_path = "java",
                             picard_path = "picard.jar"){

  #======sanity checks=======
  msg <- character()

  if(!.check_system_install("samtools")){
    msg <- c(msg, "No SAMtools install located on system path.\n")
  }
  if(!.check_system_install("bwa")){
    msg <- c(msg, "No SAMtools install located on system path.\n")
  }

  if(length(msg) != 0){
    stop(msg)
  }

  #======run================
  # check for bgzip

  cmd <- paste0("htsfile ", reference)
  res <- system(cmd, intern = TRUE)

  if(!grepl("BGZF-compressed sequence data", res)){
    if(grepl("\\.gz$", reference)){
      system(paste0("gunzip ", reference))
      system(paste0("bgzip ", gsub("\\.gz$", "", reference)))
      rm(paste0(reference, c(".amb", ".ann", ".fai", ".bwt", ".gzi", ".pac", ".sa")))
      rm(paste0(tools::file_path_sans_ext(reference), ".dict"))
    }
    else{
      system(paste0("bgzip ", reference))
      reference <- paste0(reference, ".gz")
    }
  }

  # index -- bwa
  target_files <- paste0(reference, c(".bwt", ".amb", ".ann", ".pac", ".sa"))
  if(any(!file.exists(target_files))){
    system(paste0("bwa index ", reference))
  }

  # index -- picard
  if(!file.exists(paste0(tools::file_path_sans_ext(reference, compression = TRUE), ".dict"))){
    cmd <- paste0(java_path, " -jar ", picard_path,
                  " CreateSequenceDictionary R= ", reference)
    system(cmd)
  }

  # index -- faidx
  if(!file.exists(paste0(reference, ".fai"))){
    system(paste0("samtools faidx ", reference))
  }

  return(reference)
}

#' Concatenate ".vcf" files.
#'
#' Concatenates (joins) ".vcf" files with different loci but identical individuals
#' into a single ".vcf" file.
#'
#' @param vcfs Character, vector of "vcf" files.
#' @param outfile Character,
#'   default \code{paste0(tools::file_path_sans_ext(vcfs[1], "_concat.vcf"))}.
#'   Name for final, concatenated "vcf" file.
#' @export
concat_vcfs <- function(vcfs, outfile = paste0(tools::file_path_sans_ext(vcfs[1]), "_concat.vcf")){
  #===sanity checks======
  if(!.check_system_install("bcftools")){
    stop("No bcftools install located on system path.\n")
  }

  msg <- character()

  vcfs <- normalizePath(vcfs)
  bad_vcfs <- which(!file.exists(vcfs))
  if(length(bad_vcfs) != 0){
    msg <- c(msg, paste0("Cannot locate vcfs: ", paste0(bad_vcfs, collapse = ", "), "\n"))
  }

  if(length(msg) != 0){
    stop(msg)
  }

  #========run=========
  cmd <- paste0("bcftools concat ", vcfs, " > ", outfile)
  system(cmd)

  return(outfile)
}


#' Genotype bamfiles with GATK
#'
#' Run a HaplotypeCaller, GenomicsDBI, GenotypeGVCFs, VariantFiltration
#' pipeline to genotype and hardfilter bamfiles with GATK. GATK is much more
#' suitable for higher-coverage whole-genome sequencing data than ANGSD. Can
#' have \emph{very} long run-times.
#'
#' Genotyping with GATK is a multistep, relatively complex process.
#' \code{genotype_bams_GATK} implements a pipeline to run through the major
#' steps in genotyping and filtering starting from aligned bamfiles, which can
#' be generated via \code{\link{align_denovo}} or \code{\link{align_reference}}.
#' This pipeline is \emph{much} slower than genotyping via
#' \code{\link{genotype_bams_ANGSD}}, but generally performs better for
#' high-coverage (7x+) data, particularly high-coverage whole-genome sequencing
#' data.
#'
#' For more fine-grain parallel control, the individual component functions
#' \code{prep_genome_GATK}, \code{add_RGS}, \code{run_HaplotypeCaller},
#' \code{run_GenomicsDBImport}, \code{run_GenotypeGVCFs}, and
#' \code{run_VariantFiltration} can be run piecemeal. Note that \code{add_RGs}
#' adds basic read-group information assuming each bamfile has a unique read
#' group from the same library--manually adding RGs with \code{picard}`'s
#' \code{AddOrReplaceReadGroups} utility might be preferable if more detailed
#' RG info is available, in which case the \code{add_RGs} argument should be set
#' to \code{FALSE}.
#'
#' @param bamfiles Character, vector of bamfile paths. Bamfiles must be indexed,
#'   as by default if produced using \code{\link{align_denovo}} or
#'   \code{\link{align_reference}}.
#' @param reference Character, path to a reference genome. Must be indexed
#'   with \code{bwa} and \code{picard} and zipped with \code{samtools bgzip}.
#'   If running with \code{genotype_bams_GATK}, this will happen
#'   automatically, if not, consider running \code{prep_genome_GATK}.
#' @param fastqs Character, paths to fastq files in the same order as bamfiles
#'   (if required); for \code{add_RGs}. Not required if \code{add_RGs} is
#'   \code{FALSE}.
#' @param par Numeric, default 1. Number of cores to use for parallel processing.
#' @param min_chr_size Numeric, default 0. Minimum size for chromosomes/scaffolds
#'   for genotype calling. Useful to speed up computation by excluding small
#'   unplaced scaffolds.
#' @param chunk_size Numeric or "chr", default "chr". Size of chunks by which
#'   to split up databases and genotyping, in base pairs. If "chr", chromosomes
#'   will not be split.
#' @param QD Numeric, default 2. Minimum quality-by-depth score to keep
#'   variants. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param FS Numeric, default 60. Fisher strand bias score above which variants
#'   will filtered out. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param SOR Numeric, default 3. Strand odds-ratio test score test value above
#'   which variants will be filtered out. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param MQ Numeric, default 40. Mapping quality score below which variants
#'   will be filtered out. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param MQRankSum Numeric, default -12.5. Mapping quality rank-sum test score
#'   below which variants will be filtered out. More negative values indicate
#'   that the reads with the reference alleles map better than the alternative,
#'   positive values indicate the inverse. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param ReadPosRankSum Numeric, default -8. Minimum read position rank-sum test
#'   score, below which variants will be filtered out. Measures the likelihood
#'   that alternative alleles are more likely to appear farther along in reads
#'   than the reference alleles. See
#'   \href{https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants}{here}
#'   for details.
#' @param min_genotype_quality Numeric, default 13. Phred-scaled genotype
#'   quality score below which individual genotypes will be marked as missing.
#'   The default corresponds to ~95% certainty of correct genotyping, which is
#'   usually in the approximately 7x coverage range.
#' @param platform Character, default "ILLUMINA". The platform to mark for
#'   read-groups, ignored if \code{add_RGs} is \code{FALSE} for
#'   \code{genotype_bams_GATK}.
#' @param batch_size Numeric, default 5. The number of "gvcf" files to import at
#'   a time during \code{GenomicsDBImport} operation. Higher numbers are faster
#'   but more memory intense, for the most part.
#' @param add_RGs Logical, default \code{TRUE}. If \code{TRUE}, basic read-group
#'   information will be added to bam files using the function \code{add_RGs}.
#'   If more complicated read-group information is already present in bam files,
#'   this should be set to \code{FALSE}.
#' @param java_path Character, default "java". Path to java install, by default
#'   expects java on system path.
#' @param gatk4_path Character, default "gatk.jar". Path to gatk4 ".jar" file.
#'   Used over a GATK direct system path install to allow for manual memory
#'   allocation.
#' @param picard_path Character, default "gatk.jar". Path to picard ".jar" file.
#' @param concatenate_final_vcfs Character, default \code{TRUE}. If \code{TRUE},
#'   the final "vcf" files produced by \code{GenotypeGVCFs} will be concatenated
#'   automatically using \code{\link{concat_vcfs}}.
#'
#' @export
#' @author William Hemstrom
#' @return A vector of final ".vcf"/".gvcf"/"_db"/".bam"/.etc file paths,
#'   depending on the function function.
genotype_bams_GATK <- function(bamfiles, reference, fastqs = NULL, par = 1, min_chr_size = 0, chunk_size = "chr",
                               QD = 2,
                               FS = 60,
                               SOR = 3,
                               MQ = 40,
                               MQRankSum = -12.5,
                               ReadPosRankSum = -8,
                               min_genotype_quality = 13,
                               platform = "ILLUMINA",
                               batch_size = 5,
                               add_RGs = TRUE,
                               java_path = "java", gatk4_path = "gatk.jar", picard_path = "picard.jar",
                               concatenate_final_vcfs = TRUE){

  #============sanity checks==================
  msg <- character()
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

  if(add_RGs){
    fastqs <- normalizePath(fastqs)
    bad_fastqs <- which(!file.exists(fastqs))
    if(length(bad_bams) != 0){
      msg <- c(msg, paste0("Cannot locate fastqs: ", paste0(bad_fastqs, collapse = ", "), "\n"))
    }
  }


  reference <- normalizePath(reference)
  if(!file.exists(reference)){
    msg <- c(msg, paste0("Cannot locate reference at: ", reference, "\n"))
  }

  if(min_chr_size < 0){
    msg <- c(msg, "The minimum chr size must 0 or greater.\n")
  }

  if(length(msg) > 0){
    stop(msg)
  }

  #==========setup==========
  # check ref
  reference <- prep_genome_GATK(reference, java_path, picard_path)

  # add RGs
  if(add_RGs){
    bamfiles <- add_RGS(bamfiles = bamfiles,
                        fastqs = fastqs,
                        par = par,
                        platform = platform,
                        java_path = java_path, picard_path = picard_path)
  }

  bedfiles <- make_region_beds(reference, min_chr_size = min_chr_size, chunk_size, outdir = "./")

  #=========run pipeline=====
  hapmap <- run_HaplotypeCaller(bamfiles, reference,
                                mem = mem,
                                par = par,
                                java_path = java_path, gatk4_path = gatk4_path)

  genome_dbs <- run_GenomicsDBImport(hapmap,
                                     bedfiles,
                                     mem = mem,
                                     par = par,
                                     batch_size = batch_size,
                                     java_path = java_path, gatk4_path = gatk4_path)

  raw_vcfs <- run_GenotypeGVCFs(bedfiles, reference,
                                mem = mem,
                                par = par,
                                java_path = java_path, gatk4_path = gatk4_path)

  vcfs <- run_VariantFiltration(raw_vcfs, reference,
                                mem = mem,
                                QD = QD, FS = FS, SOR = SOR, MQ = MQ, MQRankSum = MQRankSum, ReadPosRankSum = ReadPosRankSum,
                                min_genotype_quality = min_genotype_quality,
                                java_path = java_path, gatk4_path = gatk4_path)

  # concat if requested
  if(concatenate_final_vcfs & length(vcfs) > 1){
    vcfs <- concat_vcfs(vcfs)
  }

  return(vcfs)
}
