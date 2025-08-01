# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(bamfiles = params$.depends[1],
                     region = params$.depends[2]),
                params)

print(run_params)
alignR:::.rule_caller_w_slurm(alignR::run_HaplotypeCaller, run_params, "HaplotypeCaller")
cat(alignR:::.console_hline(), "\n")
