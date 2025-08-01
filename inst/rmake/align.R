# reminder, this script will be source()'d. It gets passed a named list `params` from the host function.
run_params <- c(list(RA_fastqs = params$.depends[1],
                     RB_fastqs = params$.depends[2]),
                params)

run_params$read_metadata <- run_params$read_metadata[run_params$read_metadata$sn == params$.depends[1],,drop=FALSE]
run_params$read_metadata <- run_params$read_metadata[,-1, drop = FALSE ]

print(run_params)

cat(alignR:::.console_hline(), "\n")
alignR:::.rule_caller_w_slurm(alignR::align_reference, run_params, "align")

# rename files
in_names <- paste0(rmake::replaceSuffix(params$.depends[1], ""),
                   c(".sort.flt.bam",
                     ".sort.flt.bam.bai",
                     ".sort.flt.bam.flagstat",
                     ".sort.flt.bam.stats",
                     ".markdup.bam.flagstat"))
file.rename(in_names,
            params$.target)

