\name{read2_paths_ws}
\alias{read2_paths_ws}
\title{
Read Pathways ; Map SNPs to gene-sets/pathways
}
\description{
Reads the gene-sets and uses "genome_order output"(ordered_alldata, gs_locs).\cr 
Maps the SNPs to the gene-sets
}
\usage{
read2_paths_ws(ordered_alldata = ordered_alldata, gs_locs = gs_locs, 
pathinfo = pathways)
}
\arguments{
  \item{ordered_alldata}{
Ordered Data according it the SNPs genomic location. Traits start at column 7 \cr
Return variable from:\cr 
genome_results <-genome_order(all_data=all_data)\cr
ordered_alldata <- genome_results$ordered_alldata\cr
}
  \item{gs_locs}{
Genes Annotation,Indexes and Number of observations
Return variable from:\cr
genome_results <-genome_order(all_data=all_data)
gs_locs <- genome_results$gs_locs
}
  \item{pathinfo}{
Data file with pathway id, description, and number of genes in the pathway.
PathId must correspond to the same ID of the pathways file names
 }
}
\value{
\item{pathways2}{
Pathway Id, Description,Number of Genes in the pathway, Number of genes found from the dataset,
Number of SNPs found from the dataset}
\item{per_ors}{A list of identifiers mapped to each pathway}
}

\section{Format}{
\preformatted{
Input: Ordered_alldata
name     Chromosome  Location1st GENE_ID   Symbol Orientation abpi abpilba
rs1001567          1     9194614    <NA>     <NA>        <NA> 0.96 0.89
rs1000313          1    15405489   23254 KIAA1026           + 0.93 0.57
rs1002365          1    19797248    <NA>     <NA>        <NA> 0.68 0.58
rs1002706          1    25051153    <NA>     <NA>        <NA> 0.71 0.02
rs1002487          1    26865971    6195  RPS6KA1           + 0.98 0.78

Input:gs_locs
      Symbol   Chromosome Location    Gene_ID Start_Indx Observations
 [1,] "ACYP2"  "2"        "54399633"  "98"    "35"       "1"         
 [2,] "AMPD3"  "11"       "10514707"  "272"   "898"      "1"         
 [3,] "ANK2"   "4"        "113830885" "287"   "479"      "4"         

Input:pathinfo
      PathID     PathName                        GenesCount
 [1,] "hsa00010" "Glycolysis / Gluconeogenesis"  "  66"    
 [2,] "hsa00020" "Citrate cycle (TCA cycle)"     "  31"    
 [3,] "hsa00030" "Pentose phosphate pathway"     "  27"

Output:pathways2 
ID         Name                           GenesInPath GenesFound SNPsInPath
"hsa00010" "Glycolysis / Gluconeogenesis" "  66"      "1"        "1"       
"hsa00020" "Citrate cycle (TCA cycle)"    "  31"      "0"        "0"       
"hsa00030" "Pentose phosphate pathway"    "  27"      "1"        "1"       
}
}
\seealso{
genome_permutation
}
\examples{
##DEMO / WORKSPACE#############################
data(demo,SNPs2Genes_ALL,SNPs2Genes_Distance0)
all_data <- read_pvals_ws(data_name=demo,snps_ann=SNPs2Genes_Distance0,
snps_ann_all=SNPs2Genes_ALL)
genome_results <-genome_order(all_data=all_data)
ordered_alldata <- genome_results$ordered_alldata
gs_locs <- genome_results$gs_locs
data(pathways)
data(hsa00010,hsa00020,hsa00030,hsa00040,hsa04144)
data(hsa04510,hsa04740,hsa04810,hsa05016,hsa05200)	

### Read 2 Paths from workspace
paths_res <- read2_paths_ws(ordered_alldata=ordered_alldata,
gs_locs=gs_locs,pathinfo=pathways)
#####################################
}
\keyword{ pathways }


