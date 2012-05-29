\name{demo}
\alias{demo}
\docType{data}
\title{
GWAS p_values demo data
}
\description{
GWAS p_values (tab delimited file). First Column must contain the SNP ids and the column name = "name"
}
\usage{data(demo)}
\format{
  A data frame with SNPs identifiers and gwas p-values of association
  \describe{
    \item{\code{name}}{a character vector}
    \item{\code{abpi}}{a numeric vector}
    \item{\code{abpilba}}{a numeric vector}
    \item{\code{abpildfa}}{a numeric vector}
    \item{\code{abpilpta}}{a numeric vector}
    \item{\code{abpirba}}{a numeric vector}
    \item{\code{abpirdfa}}{a numeric vector}
    \item{\code{abpirpta}}{a numeric vector}
    \item{\code{alb}}{a numeric vector}
    \item{\code{avdbp}}{a numeric vector}
  }
  \preformatted{
	name       abpi      abpilba    abpildfa  abpilpta  abpirba   abpirdfa
	rs10000010 0.9122360 0.30088096 0.2332038 0.5193068 0.1255104 0.07253145 
	rs10000023 0.8642906 0.52064064 0.9243443 0.7177759 0.9512171 0.81716250 
	rs10000030 0.2832705 0.99021664 0.8359339 0.9662707 0.8491221 0.50208681 

  }
}
\examples{
data(demo)
## use: input file for "read_pvals_fun" function
}
\keyword{datasets}
