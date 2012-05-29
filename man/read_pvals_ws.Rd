\name{read_pvals_ws}
\alias{read_pvals_ws}
\title{
Read GWAS p-values of association
}
\description{
Reads and Formats GWAS p-values of association for analysis
}
\usage{
read_pvals_ws(data_name = "", snps_ann = "", snps_ann_all = "")
}
\arguments{
  \item{data_name}{
GWAS p_values (tab delimited file)(SNP_IDs Trait1 Trait2 ...TraitN)
}
  \item{snps_ann}{
	SNPs Annotation at desired distance (SNPs2Genes_Distance0.txt)\cr
}
  \item{snps_ann_all}{
	SNPs to Gene Annotation (Complete set) (SNPs2Genes_ALL.txt) \cr
}
}
\value{
Dataframe:name;chromosome;Location;GeneID;Symbol;Orientation;Trait1;Trait2;TraitN
}
\section{Formats}{
\preformatted{
GWAS p_values (tab delimited file)(SNP_IDs Trait1 Trait2 ...TraitN)
name       abpi      abpilba    abpildfa  
rs10000010 0.9122360 0.30088096 0.2332038
rs10000023 0.8642906 0.52064064 0.9243443 
rs10000030 0.2832705 0.99021664 0.8359339

SNPs Annotation at desired distance (SNPs2Genes_Distance0.txt)
name        Chromosome  Location  GENE_ID   Symbol    Orientation
rs1000313   1           15405489  23254     KIAA1026  +
rs1000533   1           168282491 9095      TBX19     +
rs1000731   1           231963491 27185     DISC1     +

SNPs to Gene Annotation (Complete set) (SNPs2Genes_ALL.txt)
name          Chromosome_All Location_All
rs10000010    4              21618674
rs10000023    4              95733906
rs10000030    4              103374154
rs1000007     2              237752054

Output:
name       Chromosome  Location GENE_ID Symbol Orientation      abpi  
rs10000010          4  21618674   80333 KCNIP4           - 0.9122360
rs10000023          4  95733906     658 BMPR1B           + 0.8642906
rs10000030          4 103374154      NA   <NA>        <NA> 0.2832705 
}
}

\seealso{
genome_order
}
\examples{
## DEMO // WORKSPACE
data(demo,SNPs2Genes_ALL,SNPs2Genes_Distance0)
all_data <- read_pvals_ws(data_name=demo,snps_ann=SNPs2Genes_Distance0,
snps_ann_all=SNPs2Genes_ALL)
}
\keyword{ gwas }
\keyword{ pvalues }
