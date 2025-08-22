#' Prepare a conda environment with `alignR`'s dependencies
#'
#' Create a conda environment which has all of `alignR`'s main dependencies
#' installed save `alignR` itself which can be installed by activating the
#' env and installing `alignR` as usual.
#'
#' @param name character, default "alignR". The name of the conda env to create.
#'
#' @export
#' @author William Hemstrom
#' @examples
#' \dontrun{
#' # prepare environment
#' install_dependencies_conda("alignR")
#'
#' # in terminal:
#' conda activate alignR
#' R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
#' R -e "remotes::install_github('hemstrow/alignR')"
#'
#' }
install_dependencies_conda <- function(name = "alignR"){
  # if(!.check_system_install("conda")){
  #   # install conda?
  #   stop("Conda install not detected.\n")
  # }

  if(reticulate::condaenv_exists(name)){
    cat("Designated env already exists. Install dependencies to that env? y/n\n")
    resp <- readLines()
    while(resp != "y"){
      resp <- tolower(resp)
      if(resp == "yes") resp <- "y"
      if(resp == "no") resp <- "n"

      if(resp == "n"){
        stop("Designated env already exists.\n")
      }
      cat("Designated env already exists. Install dependencies to that env? y/n\n")
      resp <- readLines()
    }
  }

  reticulate::conda_create(name)
  reticulate::conda_install(name, c("r-base", "perl"))
  reticulate::conda_install(name, c("bwa",
                                    "vcftools",
                                    "bcftools",
                                    "samtools",
                                    "stacks",
                                    "gatk4",
                                    "angsd",
                                    "picard",
                                    "fastp"), channel = "bioconda")
  # cmd <- paste0("conda run -n ", name, " ", r"(R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org');remotes::install_github('hemstrow/alignR');")")

  message("Conda environment: '", name, "' prepared with `alignR` dependencies. To load, type `conda activate ", name, "`, run R, and then install `alignR`.\n")
  message("Conda env prepared with dependencies and can be activated with:\n\tconda activate ", name, "\nTo finish env, run this to install alignR:\n\t",
          paste0("conda run -n ", name, " ",
                 r"(R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org');remotes::install_github('hemstrow/alignR');")"))
  return(name)
}
