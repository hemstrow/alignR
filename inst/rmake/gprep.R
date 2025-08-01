# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
alignR:::.rule_caller_w_slurm(function(ref) {
  # index
  cmd <- paste0("bwa index ", ref)
  system(cmd)
  cmd <- paste0("samtools faidx ", ref)
  system(cmd)
  ref_bn <- gsub("\\.gz$", "", ref)
  ref_bn <- rmake::replaceSuffix(ref_bn, "")
  cmd <- paste0("picard CreateSequenceDictionary R= ", ref, " ",
                "O= ", paste0(ref_bn, ".dict"))
  system(cmd)
},
list(ref = params$.depends[1]), "genome_index")
