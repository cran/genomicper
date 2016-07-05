\name{read2_paths}
\alias{read2_paths}
\title{
Read to SNPs to sets; Map SNPs to gene-sets/pathways
}
\description{
Reads the sets/pathways, map the SNPs and genes to the gene-sets/pathways
read2_paths uses the "genome_order" output(ordered_alldata, gs_locs) to assign genomic location indexes
to each element in the gene-set. The permutation method must be defined (i.e. level = "snp" OR level = "gene").

}
\usage{
read2_paths(ordered_alldata="",gs_locs="",sets_from="workspace",
sets_prefix="hsa",level="snp",envir="")
}
\arguments{
  \item{ordered_alldata}{
 Ordered data according to the SNPs genomic location. Traits start at column 7\cr
 Return variable from:\cr
 	genome_results <-genome_order(all_data=all_data)\cr
 	ordered_alldata <- genome_results$ordered_alldata\cr
}
  \item{gs_locs}{
 Gene annotation,indexes and number of observations\cr
 Return variable from genome_order():\cr
 	genome_results <-genome_order(all_data=all_data)\cr
 	gs_locs <- genome_results$gs_locs\cr
}
  \item{sets_from}{
 Location of the gene-sets. Default set to "workspace" \cr
 sets_from="workspace" OR sets_from="directory"\cr
 "directory", only will search for information in the working directory.
}
  \item{sets_prefix}{
 Prefix of the gene-set variables or files.\cr
 Default set to sets_prefix= "hsa" 
 e.g. Variables "hsa00010","hsa00020". OR files "hsaXXXXX.txt"\cr
 each variable/file contains the list of gene identifiers part of that pathway
}
  \item{level}{
 The level at which the permutations will be performed. Assigns the indexes according to snps or genes\cr
 Default value "snp" 
 level values = "snp" OR "gene"
}
  \item{envir}{
 R environment where pathway data is stored. e.g(envir=.GlobalEnv, envir=gper.env)
}
}
\value{
\item{pathways}{
 Pathway Id, Description,Number of Genes in the pathway, Number of genes found in the dataset,
 Number of SNPs found in the dataset}
\item{per_ors}{
 A list of identifiers mapped to each pathway
 }
}
\section{Format}{
\preformatted{
Input: Ordered_alldata
name     Chromosome  Location GENE_ID   Symbol Orientation abpi abpilba
rs1001567          1 9194614     <NA>     <NA>        <NA> 0.96 0.89
rs1000313          1 15405489   23254 KIAA1026           + 0.93 0.57
rs1002365          1 19797248    <NA>     <NA>        <NA> 0.68 0.58
rs1002706          1 25051153    <NA>     <NA>        <NA> 0.71 0.02
rs1002487          1 26865971    6195  RPS6KA1           + 0.98 0.78

Input:gs_locs
      Symbol   Chromosome Location    Gene_ID Start_Indx Observations
 [1,] "ACYP2"  "2"        "54399633"  "98"    "35"       "1"         
 [2,] "AMPD3"  "11"       "10514707"  "272"   "898"      "1"         
 [3,] "ANK2"   "4"        "113830885" "287"   "479"      "4"
 
Input:pathway example
hsa04720
 [1]  10411    107  11261    114   1387 163688   ....
          
Output:pathways
ID         Name                           GenesInPath GenesFound SNPsInPath
"hsa00010" "Glycolysis / Gluconeogenesis" "  66"      "1"        "1"       
"hsa00020" "Citrate cycle (TCA cycle)"    "  31"      "0"        "0"       
"hsa00030" "Pentose phosphate pathway"    "  27"      "1"        "1"       
}
}
\seealso{
\code{\link{genes_permutation}}
\code{\link{snps_permutation}}
\code{\link{genome_order}}
}
\examples{
## DEMO - SNP Level data stored in workspace #######################
# library(genomicper)
data(demo,SNPsAnnotation)
all_data <- read_pvals(data_name=demo,snps_ann=SNPsAnnotation)
genome_results <-genome_order(all_data=all_data)
ordered_alldata <- genome_results$ordered_alldata
gs_locs <- genome_results$gs_locs
data(hsa00100,hsa00120,hsa00130,hsa00140,hsa00190,hsa02010)

paths_res <- read2_paths(ordered_alldata=ordered_alldata,
gs_locs=gs_locs,sets_from="workspace",sets_prefix="hsa",
level="snp",envir=.GlobalEnv)		
pers_ids <- paths_res$per_ors
pathways<- paths_res$pathways
####################################################################
}
\keyword{ gene-set }
\keyword{ pathways }
\keyword{ annotation }
