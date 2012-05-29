\name{genomicper-package}
\alias{genomicper-package}
\alias{genomicper}
\docType{package}
\title{
Circular Genomic Permutations
}
\description{
Circular genomic permutation uses GWAS results to establish the significance of pathway associations whilst accounting for genomic structure. All SNPs in the GWAS are placed in a 'circular genome' according to their location. Then the complete set of SNP association p-values are permuted by rotation with respect to the SNPs' genomic locations.}
\details{
\tabular{ll}{
Package: \tab genomicper\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2012-05-29\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Claudia Cabrera, Pau Navarro, Chris Haley \cr
Maintainer: Claudia Cabrera <claudia.cabrera@ed.ac.uk>
}
\references{
Uncovering Networks from Genome-Wide Association Studies via\cr 
Circular Genomic Permutation\cr
Claudia P. Cabrera*, Pau Navarro*, Jennifer E. Huffman, Alan F. Wright, Caroline Hayward,Harry Campbell, James F. Wilson, Igor Rudan, Nicholas D. Hastie, Veronique Vitart, Chris S. Haley*
}
\keyword{ package }
\examples{

### DEMO / WORKSPACE ############################################################
# library(genomicper)
# load data (SNPs GWAS, SNPs to Gene Annotations, and SNPs genomic locations)
data(demo,SNPs2Genes_Distance0,SNPs2Genes_ALL)

# format data
all_data<-read_pvals_ws(data_name=demo,snps_ann=SNPs2Genes_Distance0,
snps_ann_all=SNPs2Genes_ALL)

# order according to their genomic region
genome_results <-genome_order(all_data=all_data)

# return variables from ordering dataset
ordered_alldata <- genome_results$ordered_alldata
gs_locs <- genome_results$gs_locs

# load pathway information
data(pathways)

# load pathways, each variable contains a list of gene identifiers, 
# pathway identifiers should be the SAME identifier 
# of the pathway information file.
data(hsa00010,hsa00020,hsa00030,hsa00040,hsa04144,hsa04510,hsa04740,hsa04810)
data(hsa05016,hsa05200)

# Map SNPs to gene-sets/pathways	
paths_res <- read2_paths_ws(ordered_alldata=ordered_alldata,gs_locs=gs_locs,
pathinfo=pathways)

# Perform permutations FOLLOWING FUNCTION SAVES PERMUTATION RESULTS TO WORKSPACE 
results <- genome_permutation_ws(ordered_alldata=ordered_alldata,
pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits=c(7,10),
nper=5,threshold=0.05)
# Summary Results 
hypers<-get_results_ws()

###  ALTERNATIVE function to perform permutations:
###  save results to WORKING DIRECTORY:
# results <- genome_permutation(ordered_alldata=ordered_alldata,
# pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits=c(7,10),
# nper=5,threshold=0.05)
# hypers<-get_results(res_pattern="Permus")
####################################################################################
}
