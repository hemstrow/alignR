# start-up checks
.onAttach <- function(libname, pkgname){

  packageStartupMessage("alignR requires a few UNIX specific dependencies.\n")
  sys <- Sys.info()["sysname"]
  if(sys == "Windows"){
    packageStartupMessage("Some of these can be tricky to install on windows. We suggest using a Linux virtual machine (check out Windows Subsystem for Linux, for example). A walkthough to getting alignR up and running on Windows is availablein the 'alignR_setup' vignette.\n")
  }
  else{
    packageStartupMessage("The install_dependencies_conda() function can create a conda environment for you with all installed dependencies.\n")
  }

  # try conda
  conda_path <- try(reticulate::conda_binary(), silent = TRUE)
  Sys.setenv(conda_install = FALSE)
  if(!is(conda_path, "try-error")){
    conda_check <- system(paste0(conda_path, " --version"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    if(conda_check == 0){
      paste0("conda is good-to-go... all other dependencies can be installed into an environment with `make_conda_env()`\n")
      Sys.setenv(conda_install = TRUE)
    }
  }
  if(!.check_system_install("conda")){
    warning(paste0("No conda installation detected! Install conda to allow for easy alignR dependency installation.\n"))
  }

  # try perl
  perl_check <- suppressWarnings(system("perl -v", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(perl_check == 127){
    Sys.setenv(perl_install = FALSE)
    warning(paste0("No perl installation detected! Some functions will fail:\n\t",
                   paste0(unlist(.dependency_function_match("perl")), collapse = "\n\t"),
                   "\n"))
  }
  else{
    packageStartupMessage("perl is good-to-go!\n")
    Sys.setenv(perl_install = TRUE)
  }

  # try bash
  bash_check <- suppressWarnings(system("bash --version", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(bash_check == 127){
    Sys.setenv(bash_install = FALSE)
    warning("No bash installation detected! All functions will fail.\n")
  }
  else{
    packageStartupMessage("bash is good-to-go!\n")
    Sys.setenv(bash_install = TRUE)
  }

  # try samtools
  samtools_check <- suppressWarnings(system("samtools", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(samtools_check == 127){
    Sys.setenv(samtools_install = FALSE)
    packageStartupMessage(paste0("No SAMtools installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("samtools")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("SAMtools is good-to-go!\n")
    Sys.setenv(samtools_install = TRUE)
  }

  # try bcftools
  bcftools_check <- suppressWarnings(system("bcftools", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(bcftools_check == 127){
    Sys.setenv(bcftools_install = FALSE)
    packageStartupMessage(paste0("No bcftools installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("bcftools")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("bcftools is good-to-go!\n")
    Sys.setenv(bcftools_install = TRUE)
  }


  # try vcftools
  vcftools_check <- suppressWarnings(system("vcftools", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(vcftools_check == 127){
    Sys.setenv(vcftools_install = FALSE)
    packageStartupMessage(paste0("No vcftools installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("vcftools")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("vcftools is good-to-go!\n")
    Sys.setenv(vcftools_install = TRUE)
  }

  # try angsd
  angsd_check <- suppressWarnings(system("angsd", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(angsd_check == 127){
    Sys.setenv(angsd_install = FALSE)
    packageStartupMessage(paste0("No angsd installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("angsd")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("angsd is good-to-go!\n")
    Sys.setenv(angsd_install = TRUE)
  }


  # try bwa
  bwa_check <- suppressWarnings(system("bwa", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(bwa_check == 127){
    Sys.setenv(bwa_install = FALSE)
    packageStartupMessage(paste0("No bwa installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("bwa")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("bwa is good-to-go!\n")
    Sys.setenv(bwa_install = TRUE)
  }

  # try ngsParalog
  ngsParalog_check <- suppressWarnings(system("ngsParalog", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(ngsParalog_check == 127){
    Sys.setenv(ngsParalog_install = FALSE)
    packageStartupMessage(paste0("No ngsParalog installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("ngsParalog")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("ngsParalog is good-to-go!\n")
    Sys.setenv(ngsParalog_install = TRUE)
  }

  # try stacks
  stacks_check <- suppressWarnings(system("stacks", ignore.stdout = TRUE, ignore.stderr = TRUE))
  stacks_check2 <- suppressWarnings(system("cstacks", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(any(c(stacks_check, stacks_check2) %in% c(0, 1))){ # because an error code of 1 is actually OK too. 127 means couldn't find, a no-go
    Sys.setenv(stacks_install = ifelse(stacks_check %in% c(0, 1), "general", "specific"))
    packageStartupMessage("STACKS is good-to-go!\n")
  }
  else{
    Sys.setenv(stacks_install = FALSE)
    packageStartupMessage(paste0("No STACKS installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("stacks")), collapse = "\n\t"),
                                 "\n"))
  }

  # try fastp
  fastp_check <- suppressWarnings(system("fastp", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(fastp_check == 127){
    Sys.setenv(fastp_install = FALSE)
    packageStartupMessage(paste0("No fastp installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("fastp")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("fastp is good-to-go!\n")
    Sys.setenv(fastp_install = TRUE)
  }

  # try picard
  picard_check <- suppressWarnings(system("picard", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(picard_check == 127){
    Sys.setenv(picard_install = FALSE)
    packageStartupMessage(paste0("No picard installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("picard")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    packageStartupMessage("picard is good-to-go!\n")
    Sys.setenv(picard_install = TRUE)
  }

  # try gatk4
  gatk4_check <- suppressWarnings(system("gatk --version", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(gatk4_check == 127){
    Sys.setenv(gatk4_install = FALSE)
    packageStartupMessage(paste0("No gatk4 installation detected! Some functions will fail:\n\t",
                                 paste0(unlist(.dependency_function_match("gatk")), collapse = "\n\t"),
                                 "\n"))
  }
  else{
    gatk4_check <- suppressWarnings(system("gatk --version", intern = TRUE, ignore.stderr = TRUE))
    version <- gatk4_check[grep("GATK", gatk4_check)]
    version <- gsub(".+v", "", version)
    version <- as.numeric_version(version)
    if(version >= as.numeric_version("4.0.0.0")){
      packageStartupMessage("gatk4 is good-to-go!\n")
      Sys.setenv(gatk4_install = TRUE)
    }
    else{
      Sys.setenv(fastp_install = FALSE)
      packageStartupMessage(paste0("gatk installation detected, but older than 4.0.0.0! Some functions will fail:\n\t",
                                   paste0(unlist(.dependency_function_match("gatk")), collapse = "\n\t"),
                                   "\n"))
    }
  }
}
