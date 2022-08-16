#' Run the alignR GUI
#'
#' Runs \code{alignR}'s Graphical User Interface (GUI) \code{\link[shiny]{shiny}}
#' App to help simplify demultiplexing, alignment, and genotyping using the full
#' suite of tools available in the package.
#' 
#' @author William Hemstrom
#' @export
alignR_gui <- function(){
  if(!.check_system_install("perl")){
    msg <- c(msg, "No perl install located on system path.\n")
  }
  if(!.check_system_install("bash")){
    msg <- c(msg, "No bash install located on system path.\n")
  }
  
  missing_depends <- character()
  if(!.check_system_install("samtools")){
    missing_depends <- c(missing_depends, "samtools")
  }
  if(!.check_system_install("angsd")){
    missing_depends <- c(missing_depends, "angsd")
  }
  if(!.check_system_install("ngsParalog")){
    missing_depends <- c(missing_depends, "ngsParalog")
  }
  if(isFALSE(.check_system_install("stacks"))){
    missing_depends <- c(missing_depends, "stacks")
  }
  if(!.check_system_install("bwa")){
    missing_depends <- c(missing_depends, "bwa")
  }
  if(!.check_system_install("bcftools")){
    missing_depends <- c(missing_depends, "bcftools")
  }
  
  if(length(missing_depends) != 0){
    message("Note: some dependencies are missing. Dependant features will be disabled. See package start-up message for details.")
  }
  
  gui_path <- system.file("shiny", "ui_app", package = "alignR")
  if(gui_path == ""){stop("Could not locate GUI app. Package may not be properly installed.")}
  shiny::runApp(gui_path)
}