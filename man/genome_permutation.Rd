\name{genome_permutation}
\alias{genome_permutation}
\title{
Genomic Permutations
}
\description{
Performs genomic permutations (Saves Permutation results to the working directory)
}
\usage{
genome_permutation(ordered_alldata = "",pers_ids = paths_res$per_ors,
pathways2=paths_res$pathways2,ntraits = "", nper = 100, threshold = 0.05)
}
\arguments{
  \item{ordered_alldata}{
Return variable from "genome_order". Ordered genome and trait p-values
}
  \item{pers_ids}{
Return variable "per_ors" from "read2_paths"
}
  \item{pathways2}{
Return variable "pathways2" from "read2_paths"
}
  \item{ntraits}{
Traits INDEX to be analysed. Index according to "ordered_alldata". \cr
Trait Columns index must start at 7. Example: ntraits=c(7:9), ntraits=7
}
  \item{nper}{
Number of permutations.Example: nper=100
}
  \item{threshold}{
Threshold to be set by the hypergeometric test. threshold=0.05
}
}
\value{
Returns "Permus_trait" variables (the permutation datasets to the working directory). 
}
\references{
HyperBigUniverse_v2. Imports phyper (from stats)
}
\seealso{
get_results
}
\examples{
################DEMO/WORKSPACE###################################
#library(genomicper)
data(demo,SNPs2Genes_ALL,SNPs2Genes_Distance0)
all_data <- read_pvals_ws(data_name=demo,snps_ann=SNPs2Genes_Distance0,
snps_ann_all=SNPs2Genes_ALL)
genome_results <-genome_order(all_data=all_data)
ordered_alldata <- genome_results$ordered_alldata
gs_locs <- genome_results$gs_locs
data(pathways)
data(hsa00010,hsa00020,hsa00030,hsa00040)
data(hsa04144,hsa04510,hsa04740,hsa04810,hsa05016,hsa05200)	
paths_res <- read2_paths_ws(ordered_alldata=ordered_alldata,
gs_locs=gs_locs,pathinfo=pathways)
### function to perform genomic permutations and 
### save the results of each permutation to the working directory
### results <- genome_permutation(ordered_alldata=ordered_alldata,
### pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits=c(7,10),
### nper=5,threshold=0.05)
}
\keyword{ permutation }
\keyword{ genome }