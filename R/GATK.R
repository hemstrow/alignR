run_HaplotypeCaller <- function(bamfiles, reference, mem, par = 1, java_path = "java", gatk4_path = "gatk"){
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
    for(i in 1:length(chunks[[q]])){
      # run Haplotype Caller

      ttempdir <- file.path(temp_dir, paste0("_HC_q", q))
      dir.create(ttempdir)
      cmd <- paste0("bash run_HaplotypeCaller.sh ",
                    chunks[[q]][i], " ",
                    reference, " ",
                    ttempdir, " ",
                    mem, " ",
                    java_path, " ",
                    gatk4_path, " ")

      system(cmd)
      unlink(ttempdir, recursive = TRUE)
    }
  }

  # release cores and clean up
  parallel::stopCluster(cl)

  # save a hapmap
  sample_map <- data.frame(samp = basename(tools::file_path_sans_ext(bamfiles)), file = paste0(bamfiles, ".hapcalls.gvcf.gz"))
  data.table::fwrite(sample_map, paste0(bamfiles[1], "_hapmap.txt"), sep = "\t", col.names = F, row.names = F)
  cat("hapmap saved to: ", paste0(bamfiles[1], "_hapmap.txt"))

  # return info
  return(sample_map)
}

add_RGS <- function(bamfiles, fastqs, par = 1, platform = "ILLUMINA", java_path = "java", picard_path = "picard"){
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
}

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


  bed <- data.frame(gsub(" .+", "", d), 0, as.numeric(dl))

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

run_GenomicsDBImport <- function(hapmap, bedfiles, mem, par = 1, batch_size = 5, java_path = "java", gatk4_path = "gatk"){
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
  browser()
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      ttempdir <- file.path(temp_dir, paste0("_GDBImport_q", q))
      dir.create(ttempdir)
      cmd <- paste0("bash ", script, " ",
                    hapmap, " ",
                    chunks[[q]][i], " ",
                    mem, " ",
                    ttempdir, " ",
                    par, " ",
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
  bedfiles <- basename(bedfiles)
  return(paste0(bedfiles, "_db"))
}

run_GenotypeGVCFs <- function(bedfiles, reference, mem, par = 1, java_path = "java", gatk4_path = "gatk"){
  #============sanity checks==================
  msg <- character()

  bedfiles <- normalizePath(bedfiles)
  bad_beds <- which(!file.exists(bedfiles))
  if(length(bad_beds) != 0){
    msg <- c(msg, paste0("Cannot locate bedfiles: ", paste0(bad_beds, collapse = ", "), "\n"))
  }

  # check that dbs exist
  dbs <- basename(bedfiles)
  dbs <- paste0(dbs, "_db")
  bad_dbs <- which(!dir.exists(dbs))
  if(length(bad_dbs) != 0){
    msg <- c(msg, paste0("Cannot locate GenomicsDBs: ", paste0(bad_beds, collapse = ", "), "\nDBs are expected to be named matching .bed files, ending in _db"))
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

  # run
  output <- foreach::foreach(q = 1:par, .inorder = FALSE, .errorhandling = "pass",
                             .packages = "alignR"
  ) %dopar% {
    for(i in 1:length(chunks[[q]])){
      # run Haplotype Caller

      ttempdir <- file.path(temp_dir, paste0("_GenotypeGVCFs_q", q))
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

  return(paste0(dirname(bedfiles), "/raw_", basename(tools::file_path_sans_ext(befiles)), ".vcf"))
}

run_VariantFiltration <- function(vcfs, reference, mem,
                                  QD = 2,
                                  FS = 60,
                                  SOR = 3,
                                  MQ = 40,
                                  MQRankSum = -12.5,
                                  ReadPosRankSum = -8,
                                  min_genotype_quality = 13,
                                  java_path = "java", gatk4_path = "gatk"){
  #============sanity checks==================
  msg <- character()

  vcfs <- normalizePath(vcfs)
  bad_vcfs <- which(!file.exists(vcfs))
  if(length(bad_beds) != 0){
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
      cmd <- paste0("bash run_VariantFiltration.sh ",
                    tools::file_path_sans_ext(chunks[[q]][i]), " ",
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

  return(paste0(tools::file_path_sans_ext(vcfs), ".recode.vcf"))
}

prep_genome_GATK <- function(reference, java_path = "java", picard_path = "picard"){

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
  if(!file.exists(paste0(tools::file_path_sans_ext(reference), ".dict"))){
    cmd <- paste0(java_path, " -jar ", picard_path,
                  " CreateSequenceDictionary R= ", reference,
                  " O= ", paste0(tools::file_path_sans_ext(reference), ".dict"))
    system(cmd)
  }

  # index -- faidx
  if(!file.exists(paste0(reference, ".fai"))){
    system(paste0("samtools faidx ", reference))
  }

  return(reference)
}

concat_vcfs <- function(vcfs, outfile = paste0(tools::file_path_sans_ext(vcfs[1], "_concat.vcf"))){
  #===sanity checks======
  if(!.check_system_install("bcftools")){
    stop("No bcftools install located on system path.\n")
  }

  msg <- character()

  vcfs <- normalizePath(vcfs)
  bad_vcfs <- which(!file.exists(vcfs))
  if(length(bad_bams) != 0){
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

genotype_bams_GATK <- function(bamfiles, reference, fastqs, par, min_chr_size = 0, chunk_size = "chr",
                               QD = 2,
                               FS = 60,
                               SOR = 3,
                               MQ = 40,
                               MQRankSum = -12.5,
                               ReadPosRankSum = -8,
                               min_genotype_quality = 13,
                               platform = "ILLUMINA",
                               batch_size = 5,
                               java_path = "java", gatk4_path = "gatk", picard_path = "picard",
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

  fastqs <- normalizePath(fastqs)
  bad_fastqs <- which(!file.exists(fastqs))
  if(length(bad_bams) != 0){
    msg <- c(msg, paste0("Cannot locate fastqs: ", paste0(bad_fastqs, collapse = ", "), "\n"))
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
  bamfiles <- add_RGS(bamfiles = bamfiles,
                      fastqs = fastqs,
                      par = par,
                      platform = platform,
                      java_path = java_path, picard_path = picard_path)

  bedfiles <- make_region_beds(reference, min_chr_size = min_chr_size, chunk_size, outdir = "./bedfiles")

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
  if(concatenate_final_vcfs){
    vcfs <- concat_vcfs(vcfs)
  }

  return(vcfs)
}
