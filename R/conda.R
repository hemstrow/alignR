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

  reticulate::conda_create(name)
  reticulate::conda_install(name, c("r-base", "perl"))
  reticulate::conda_install(name, c("bwa",
                                    "vcftools",
                                    "bcftools",
                                    "samtools",
                                    "stacks",
                                    "gatk4",
                                    "picard"), channel = "bioconda")

  message("Conda environment: '", name, "' prepared with `alignR` dependencies. To load, type `conda activate ", name, "`, run R, and then install `alignR`.\n")
  message("In shell, run:\nconda activate ", name, "\n",
          r"(R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org');remotes::install_github('hemstrow/alignR');")")
  return(name)
}
