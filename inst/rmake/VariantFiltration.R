# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(vcfs = params$.depends[1],
                     reference = params$.depends[2]),
                params)

print(run_params)
alignR:::.rule_caller_w_slurm(alignR::run_VariantFiltration, run_params, "HardFilter")

# rename
cat(alignR:::.console_hline(), "\n")
