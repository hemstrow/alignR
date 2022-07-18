.fetch_a_script <- function(script_name, dir) system.file(dir, script_name, "alignR")


.check_system_install <- function(code) as.logical(Sys.getenv(paste0(code, "_install")))

