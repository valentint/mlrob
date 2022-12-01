\name{disc.mentr}
\alias{disc.mentr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Discretization using the minimum entropy criterion}
\description{
  This function discretizes the continuous attributes of a data frame 
using the minumum entropy criterion along with the minimum description length 
as stopping rule.
}
\usage{
disc.mentr(data, varcon,out=c("symb","num"))
}

\arguments{
  \item{data}{ The name of the dataset to be discretized}
  \item{varcon}{ A vector containing the indices of the columms 
to be discretized}
\item{out}{To get the data discretized in numerical form enter "num".
To get the data discretized in interval form enter "symb"}
}
\value{
Returns a matrix containing only discretized features.}
\references{ 
Dougherty, J., Kohavi, R., and Sahami, M. (1995). Supervised and unsupervised discretization 
of continuous features. ML-95.
}
\author{Luis Daza}
\seealso{ \link{disc.1r}, \link{disc.ew},\link{disc.ef},\link{chiMerge} }
\examples{
#----Discretization using the entropy with Minimum Description Length.
data(bupa)
bupa.disc=disc.mentr(bupa,1:6,out="num")
}
\keyword{Discretization}
