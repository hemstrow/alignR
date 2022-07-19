# start-up checks
.onAttach <- function(libname, pkgname){

  # try perl
  perl_check <- suppressWarnings(system("perl -v", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(perl_check == 127){
    Sys.setenv(perl_install = FALSE)
    warning("No perl installation detected! Some functions will fail.\n")
  }
  else{
    packageStartupMessage("perl is good-to-go!\n")
    Sys.setenv(perl_install = TRUE)
  }

  # try bash
  bash_check <- suppressWarnings(system("bash --version", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(bash_check == 127){
    Sys.setenv(bash_install = FALSE)
    warning("No bash installation detected! Some functions will fail.\n")
  }
  else{
    packageStartupMessage("bash is good-to-go!\n")
    Sys.setenv(bash_install = TRUE)
  }

  # try samtools
  samtools_check <- suppressWarnings(system("samtools", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(samtools_check == 127){
    Sys.setenv(samtools_install = FALSE)
    packageStartupMessage("No samtools installation detected! Some functions will fail.\n")
  }
  else{
    packageStartupMessage("Samtools is good-to-go!\n")
    Sys.setenv(samtools_install = TRUE)
  }

  # try angsd
  angsd_check <- suppressWarnings(system("angsd", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(angsd_check == 127){
    Sys.setenv(angsd_install = FALSE)
    packageStartupMessage("No angsd installation detected! Some functions will fail.\n")
  }
  else{
    packageStartupMessage("angsd is good-to-go!\n")
    Sys.setenv(angsd_install = TRUE)
  }


  # try bwa
  bwa_check <- suppressWarnings(system("bwa", ignore.stdout = TRUE, ignore.stderr = TRUE))
  if(bwa_check == 127){
    Sys.setenv(bwa_install = FALSE)
    packageStartupMessage("No bwa installation detected! Some functions will fail.\n")
  }
  else{
    packageStartupMessage("bwa is good-to-go!\n")
    Sys.setenv(bwa_install = TRUE)
  }
}
