\name{discretevar}
\alias{discretevar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Performs Minimum Entropy discretization for a given attribute}
\description{
 This function carries out ME discretization for a given attribute of a 
dataset. It is also called from within the function discr.mentr.}
\usage{
discretevar(data, var, n, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The name of the dataset}
  \item{var}{ The column where the attribute to be discretized is located}
  \item{n}{ The number of rows of the dataset}
  \item{p}{ The number of columns of the dataset}
}
\author{ Luis Daza}
\seealso{\code{\link{disc.mentr}}}
\keyword{Discretization}
