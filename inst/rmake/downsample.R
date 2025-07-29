# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(bams = params$.depends[1]),
                params)
cat(alignR:::.console_hline(), "\n")
alignR:::.rule_caller_w_slurm(alignR::downsample_bams, run_params, "downsample")
