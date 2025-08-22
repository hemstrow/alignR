# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(db = params$.depends[2],
                     L = params$.depends[1],
                     outfile = params$.target[1],
                     reference = params$.depends[3]),
                params)

print(run_params)

genotype_f <- function(db, L, mem, batch_size, outfile, reference, max_alternate_alleles, new_qual = TRUE){
  temp_dir <- tempdir()

  cmd <- paste0('gatk --java-options "-Xmx', mem, 'g -Xms', mem, 'g" GenotypeGVCFs ',
                ' -V gendb://', db, ' ',
                '-R ', reference, ' ',
                '-L ', L, ' ',
                '-O ', outfile, ' ',
                '--tmp-dir ', temp_dir, ' ',
                ifelse(new_qual, '-new-qual ', ''),
                '--max-alternate-alleles ', max_alternate_alleles, ' ',
                '--disable-bam-index-caching ')

  system(cmd)
  unlink(temp_dir)
}

alignR:::.rule_caller_w_slurm(genotype_f, run_params, "align")

cat(alignR:::.console_hline(), "\n")
