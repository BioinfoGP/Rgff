% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rgff.R
\name{sort_gtf}
\alias{sort_gtf}
\title{Sorts a GTF file}
\usage{
sort_gtf(gtfFile, outFile, forceOverwrite = FALSE)
}
\arguments{
\item{gtfFile}{Path to the input GTF file}

\item{outFile}{Path to the output GTF file, if not provided the output will be  the input path (without extension) with the suffix sorted.gtf}

\item{forceOverwrite}{If output file exists, overwrite the existing file. (default FALSE)}
}
\value{
Path to the sorted GTF file
}
\description{
This function produces a sorted GTF file from an unsorted GTF file.
The default order is by Chromosome, Start, End (reverse) and feature (gene -> transcript -> OTHERS)
}
\keyword{internal}
