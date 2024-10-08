% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{merge_bams}
\alias{merge_bams}
\title{Merge bam files}
\usage{
merge_bams(file_list, names = NULL, par = 1)
}
\arguments{
\item{file_list}{list where each element contains a set of bam files to
merge. Each bam file must be sorted, see details.}

\item{names}{character vector, default NULL. A vector of names for each
merged output file. If not provided, files will be named "merged_x", where
x is the name of the first file in each set to merge.}

\item{par}{numeric, default 1. Number of cores to use for the merges.}
}
\value{
Generates merged bam files. File paths provided as a character
  vector returned from this function.
}
\description{
Merges bam files as a set from each element of a list using samtools.
Resulting merged files will be in the directory of the first bam file of the
first set to merge. Bamfiles must be sorted, see details.
}
\details{
Bamfiles must be sorted to merge them using samtools. Outputs from both
bamfile-generating alignR functions (\code{align_denovo} and
\code{align_reference}) are sorted. Otherwise, bamfiles can be sorted using
the \code{sort_bamfiles} utility function.
}
\author{
William Hemstrom
}
