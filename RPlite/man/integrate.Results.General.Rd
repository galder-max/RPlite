\name{integrate.Results.General}
\alias{integrate.Results.General}

\title{
integrate results of multiple RP outputs for meta-analysis
}

\description{
this function allows to perform posthoc meta-analysis by integrating
single dataset results together instead of recomputing them for the
meta-analysis.
Also useful for meta-analsis of mixture of datasets with paired
comparisons and datasets with unpaired comparisons.
}

\usage{
integrate.Results.General(listRes, vec = NA, weights = NA, nperms = NA)
}

\arguments{
  \item{listRes}{
    a list of the results as returned by the main RP function, i.e. RP.symetric
  }
\item{vec}{
  if listRes contains elements that are not to be included in the
  meta-analysis, vec can be specified and is the indices of the elements to
  integrate from listRes. This can be useful to select only a
  subset of datasets for the analysis.
}
  \item{weights}{
    weights can be given to each dataset. The default is the number of
  possible pairs of ratios, i.e. number of samples from level1
  multiplied by number of samples from level2 in each dataset.
}
  \item{nperms}{
    number of permutations to use. default is the lowest number of
    permutations used for the results in listRes. More can be given and
    will require more computing time.
}
}
\details{

}
\value{
 list of results with the null distribution, the final results in finalRes, and the
 final scores in final.scoresF
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}


\seealso{

}
\examples{

}
\keyword{ ~kwd1 }
