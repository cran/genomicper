\name{HyperBigUniverse_v2}
\alias{HyperBigUniverse_v2}
\title{
Hypergeometric Test (phyper)
}
\description{
Performs Hypergeometric test (phyper from R) 
}
\usage{
HyperBigUniverse_v2(Sig_in_Paths, uniSig, gns_in_Paths, universe = rowsf)
}
\arguments{
  \item{Sig_in_Paths}{
Number of significant genes in the pathway
}
  \item{uniSig}{
Number of significant genes in the dataset
}
  \item{gns_in_Paths}{
Number of genes in the pathway
}
  \item{universe}{
Number of genes in the dataset
}
}
\value{
Returns hypergeometric test
}
\references{
HyperBigUniverse_v2. Imports phyper (from stats)
}