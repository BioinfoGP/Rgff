% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rgff.R
\name{check_gtf}
\alias{check_gtf}
\title{Test consistency and order of a GTF file}
\usage{
check_gtf(gtfFile)
}
\arguments{
\item{gtfFile}{Path to the input GTF file}
}
\value{
A data frame with the errors detected
}
\description{
This function tests consistency and order of a GTF file.
Optionally, can fix some of the consistency errors and sort the GTF file.
}
\keyword{internal}
