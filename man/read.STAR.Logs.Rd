% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNASeq.metrics.R
\name{read.STAR.Logs}
\alias{read.STAR.Logs}
\title{Read log files from STAR alignment}
\usage{
read.STAR.Logs(starDir)
}
\arguments{
\item{starDir}{Character vector with the directory where STAR *Log.final.out files are stored}
}
\value{
Numeric dataframe with one row per file and columns:  number of total reads, uniquely mapped, 
multimapped, multmapped (too manany), mapped (too maney missmatches), 
unmapped (too short), unmapped (other) and chimeric.
}
\description{
Read log files from STAR alignment
}
\author{
Julia Perera Bel <jperera@imim.es>
}
