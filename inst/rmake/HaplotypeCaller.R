# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(bamfiles = params$.depends[1],
                     region = params$.depends[2]),
                params)

print(run_params)
alignR:::.rule_caller_w_slurm(alignR::run_HaplotypeCaller, run_params, "HaplotypeCaller")

# rename
file.rename(c(paste0(run_params$bamfiles, "-", basename(tools::file_path_sans_ext(run_params$region)), ".hapcalls.gvcf.gz"),
              paste0(run_params$bamfiles, "-", basename(tools::file_path_sans_ext(run_params$region)), ".hapcalls.gvcf.gz.tbi")),
            params$.target)


cat(alignR:::.console_hline(), "\n")
