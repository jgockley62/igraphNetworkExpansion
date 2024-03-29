% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{name_pull}
\alias{name_pull}
\title{Pull the names of a PathTrace}
\usage{
name_pull(char)
}
\arguments{
\item{char}{is a list entry from get.all.shortest.paths()}
}
\value{
names of the genes in the input trace path
}
\description{
Pulls the names from an igraph get.all.shortest.paths() object
}
\examples{
example_path <- list(
   c('319', "49", "23")
)
names(example_path[[1]]) <- c("GeneA", "GeneZ", "GeneAlpha")
name_pull(example_path[[1]])
}
