.fetch_a_script <- function(script_name, dir) system.file(dir, script_name, package = "alignR")

.check_system_install <- function(code) as.logical(Sys.getenv(paste0(code, "_install")))

.paired_length_check <- function(RA, RB = NULL){
  if(!is.null(RB)){
    return(length(RA) == length(RB))
  }
  else{
    return(TRUE)
  }
} 