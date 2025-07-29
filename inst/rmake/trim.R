# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(RA_fastqs = params$.depends[1],
                   RB_fastqs = params$.depends[2]),
                params)
cat(alignR:::.console_hline(), "\n")
alignR:::.rule_caller_w_slurm(alignR::trim_fastp, run_params, "fastp")
