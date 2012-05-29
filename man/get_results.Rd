\name{get_results}
\alias{get_results}
\title{
Circular Permutation Results
}
\description{
Creates a summary of the genomic permutations (permutations performed with "genome_permutation"). Permutation results stored at the directory.
}
\usage{
get_results(res_pattern = "Permus")
}
\arguments{
  \item{res_pattern}{
Pattern of the Permutation files. eg. res=pattern="Permus"
}
}
\value{
\item{hyper_vs_empirical}{
Data frame with Path ID, Trait, Threshold set by permutations,\cr 
Theoretical Hypergeometric p-value,\cr 
Observed (Empirical Hypergeometric p-values)
}
}
\section{Format}{
\preformatted{
     PathID Trait   Threshold     P-Value  Observed
1  hsa00010  abpi 0.040441176 0.058823529 1.0000000
2  hsa00020  abpi 0.000000000 0.000000000 0.1666667
3  hsa00030  abpi 0.040441176 0.058823529 1.0000000
4  hsa00040  abpi 0.004377397 0.166871669 1.0000000
5  hsa04144  abpi 0.116896390 0.166871669 1.0000000
}
}

\examples{
#### DEMO / WORKSPACE ##########################################################
data(demo,SNPs2Genes_ALL,SNPs2Genes_Distance0)
all_data <- read_pvals_ws(data_name=demo,
snps_ann=SNPs2Genes_Distance0,snps_ann_all=SNPs2Genes_ALL)
genome_results <-genome_order(all_data=all_data)
ordered_alldata <- genome_results$ordered_alldata
gs_locs <- genome_results$gs_locs
data(pathways)
data(hsa00010,hsa00020,hsa00030,hsa00040)
data(hsa04144,hsa04510,hsa04740,hsa04810,hsa05016,hsa05200)	
paths_res <- read2_paths_ws(ordered_alldata=ordered_alldata,gs_locs=gs_locs,
pathinfo=pathways)

##FROM WORKSPACE
results <- genome_permutation_ws(ordered_alldata=ordered_alldata,
pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits=c(7,10),
nper=5,threshold=0.05)
hypers<-get_results_ws(res_pattern="Permus")

##FROM DIRECTORY
# results <- genome_permutation(ordered_alldata=ordered_alldata,
# pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits=c(7,10),
# nper=5,threshold=0.05)
# hypers<-get_results(res_pattern="Permus")
}
\keyword{ results }
