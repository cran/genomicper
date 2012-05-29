\name{SNPs2Genes_ALL}
\alias{SNPs2Genes_ALL}
\docType{data}
\title{
SNPs2Genes (all genes) 
}
\description{
SNPs annotated to genes (update all SNP's genomic location)
}
\usage{data(SNPs2Genes_ALL)}
\format{
  A data frame with all SNP annotations (locations), following 3 variables.
  \describe{
    \item{\code{name}}{Character vector. SNP rs identifiers}
    \item{\code{Chromosome_All}}{Character vector, Chromosome}
    \item{\code{Location_All}}{Numeric vector, SNP location}
  }
  \preformatted{
	name       Chromosome_All Location_All
	rs10000010              4     21618674
	rs10000023              4     95733906
	rs10000030              4    103374154
	}
}
\examples{
data(SNPs2Genes_ALL)
}
\keyword{datasets}
