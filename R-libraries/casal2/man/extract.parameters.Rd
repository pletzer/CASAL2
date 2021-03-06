% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.parameters.R
\name{extract.parameters}
\alias{extract.parameters}
\title{Utility extract.parameters function}
\usage{
extract.parameters(file, path = "", fileEncoding = "")
}
\arguments{
\item{file}{the name of the input file containing model output to extract}

\item{path}{Optionally, the path to the file}

\item{fileEncoding}{Optional, allows the R-library to read in files that have been encoded in alternative UTF formats, see the manual for the error message that would indicate when to use this switch.}
}
\value{
Data <"data.frame"> of parameters that are from a -i format.
}
\description{
This function reads in a parameter file that would be generated using the -o syntax.
}
\author{
Craig Marsh
}
