% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gi2track.R
\name{gi2track}
\alias{gi2track}
\title{convert GInteractions to track object}
\usage{
gi2track(gi, gi2)
}
\arguments{
\item{gi}{an object of GInteractions}

\item{gi2}{an object of GInteractions}
}
\value{
an track object
}
\description{
Convert GInteractions object to track object
}
\examples{
gi <- readRDS(system.file("extdata", "nij.chr6.51120000.53200000.gi.rds", package="trackViewer"))
gi2track(gi)
}
