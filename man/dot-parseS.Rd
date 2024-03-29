% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/julia.R
\name{.parseS}
\alias{.parseS}
\title{internal functions to parse the results of a GenomnePermutaiton.Jl permtest
This version deals with the simple_P test}
\usage{
.parseS(results)
}
\arguments{
\item{results}{the julia object from GenomePermutation's permtest}
}
\value{
a list with the fields:
iterm, tested, randomised, obs, ran, eval_function, random_strategy, P, test, notes
}
\description{
internal functions to parse the results of a GenomnePermutaiton.Jl permtest
This version deals with the simple_P test
}
\note{
the notes field contains a line of information for plotting specific to each test
}
\seealso{
\code{\link{.parseZ}}
}
