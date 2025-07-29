# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
alignR:::.rule_caller_w_slurm(function(ref) {
  # index
  cmd <- paste0("bwa index ", ref)
  system(cmd)
  cmd <- paste0("samtools faidx ", ref)
  system(cmd)

  # length
  l <- data.table::fread(paste0(ref, ".fai"), select = 2)
  l <- sum(l[[1]])
  writeLines(as.character(l), paste0(ref, ".len"))
},
list(ref = params$.depends[1]), "genome_index")
