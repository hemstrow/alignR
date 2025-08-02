# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(gvcfs = params$.depends[grep(".gvcf.gz$", params$.depends)],
                     outfile = params$.target[1]),
                params)

print(run_params)

gdbi_f <- function(gvcfs, L, mem, batch_size, outfile){
  unlink(outfile)
  temp_dir <- tempdir()

  cmd <- paste0('gatk --java-options "-Xmx', mem, 'g -Xms', mem, 'g" GenomicsDBImport ',
                '--genomicsdb-workspace-path ', outfile, ' ',
                '--batch-size ', batch_size, ' ',
                '-L ', L, ' ',
                '--tmp-dir ', temp_dir, ' ')
  cmd <- paste0(cmd,
                paste0('-V ', gvcfs),
                ' ')

  system(cmd)
  unlink(temp_dir)
}

alignR:::.rule_caller_w_slurm(gdbi_f, run_params, "align")

cat(alignR:::.console_hline(), "\n")
