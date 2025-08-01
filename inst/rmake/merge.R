# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- params
cat(alignR:::.console_hline(), "\n")
alignR:::.rule_caller_w_slurm(alignR::merge_bams, run_params, "merge_bams")
