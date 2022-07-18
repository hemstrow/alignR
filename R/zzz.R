# start-up checks
.onAttach <- function(libname, pkgname){

  # try perl
  script <- .fetch_a_script("startup_check.pl", "perl")
  perl_check <- system(paste0("perl ", script))
  if(perl_check == 127){
    Sys.setenv(perl_install = FALSE)
    warning("No perl installation detected! Some functions will fail.\n")
  }
  else{
    Sys.setenv(perl_install = TRUE)
  }

  # try bash
  script <- .fetch_a_script("startup_check.sh", "shell")
  bash_check <- system(paste0("bash ", script))
  if(bash_check == 127){
    Sys.setenv(bash_install = FALSE)
    warning("No bash installation detected! Some functions will fail.\n")
  }
  else{
    Sys.setenv(bash_install = TRUE)
  }

  # try samtools
  samtools_check <- system("samtools", ignore.stdout = TRUE)
  if(samtools_check == 127){
    Sys.setenv(samtools_install = FALSE)
    warning("No samtools installation detected! Some functions will fail.\n")
  }
  else{
    cat("Samtools is good-to-go!\n")
    Sys.setenv(samtools_install = TRUE)
  }

  # try angsd
  angsd_check <- system("angsd", ignore.stdout = TRUE)
  if(angsd_check == 127){
    Sys.setenv(angsd_install = FALSE)
    warning("No angsd installation detected! Some functions will fail.\n")
  }
  else{
    cat("angsd is good-to-go!\n")
    Sys.setenv(angsd_install = TRUE)
  }


  # try bwa
  bwa_check <- system("bwa", ignore.stdout = TRUE)
  if(bwa_check == 127){
    Sys.setenv(bwa_install = FALSE)
    warning("No bwa installation detected! Some functions will fail.\n")
  }
  else{
    cat("bwa is good-to-go!\n")
    Sys.setenv(bwa_install = TRUE)
  }
}
