\name{SNPs2Genes_Distance0}
\alias{SNPs2Genes_Distance0}
\docType{data}
\title{
SNPS-Genes annotation to Distance 0 (within a gene)
}
\description{
SNPs annotated to genes. Annotation only when the SNPs fall within start and end of transcription of the genes.
}
\usage{data(SNPs2Genes_Distance0)}
\format{
  A data frame with 148710 observations on the following 6 variables.
  \describe{
    \item{\code{name}}{a character vector}
    \item{\code{Chromosome}}{a character vector}
    \item{\code{Location}}{a numeric vector}
    \item{\code{GENE_ID}}{a numeric vector}
    \item{\code{Symbol}}{a character vector}
    \item{\code{Orientation}}{a character vector}
  }
  \preformatted{
	name      Chromosome  Location  GENE_ID    Symbol     Orientation
	rs1000313          1  15405489    23254    KIAA1026             +
	rs1000533          1 168282491     9095    TBX19                +
	rs1000731          1 231963491    27185    DISC1                +
	}
}
\source{
NCBI Gene database,(http://www.ncbi.nlm.nih.gov/gene ; Build.37.1).
}
\examples{
data(SNPs2Genes_Distance0)
}
\keyword{datasets}
