\name{read_pvals}
\alias{read_pvals}
\title{
Read GWAS p-values of association and Merge with SNP annotations
}
\description{
Read GWAS p-values of association and Merge with SNP annotations for analysis
}
\usage{
read_pvals(data_name="",snps_ann="",from="workspace")
}
\arguments{
  \item{data_name}{
    GWAS p_values (tab delimited file)(SNP_IDs Trait1 Trait2 ...TraitN)
}
  \item{snps_ann}{
	SNPs Annotation (SNPsAnnotation). Genomicper uses entrez gene ids to annotate associate SNPs-to genes-pathways\cr
	The annotation MUST match your data input (coordinates and chromosome format)\cr
	Any SNP ID is valid, as long the ID is set as character\cr
	The examples below show an option on how to annotate the SNPs prior the use of genomicper\cr
}
  \item{from}{
	Datasets location. Values "workspace" OR "directory" \cr
}
}
\value{
Dataframe: 
name; chromosome; Location; GeneID; Symbol; Orientation; Trait1; TraitN
}
\section{Formats}{
\preformatted{
GWAS p_values (tab delimited file)(SNP_IDs Trait1 Trait2 ...TraitN)
name       abpi      abpilba    abpildfa  
rs10000010 0.9122360 0.30088096 0.2332038
rs10000023 0.8642906 0.52064064 0.9243443 
rs10000030 0.2832705 0.99021664 0.8359339

SNPs Annotation (SNPsAnnotation)
name        Chromosome  Location  GENE_ID   Symbol    Orientation
rs1000313   1           15405489  23254     KIAA1026  +
rs1000533   1           168282491 9095      TBX19     +
rs1000731   1           231963491 27185     DISC1     +

Output:
name       Chromosome  Location GENE_ID Symbol Orientation      abpi  
rs10000010          4  21618674   80333 KCNIP4           - 0.9122360
rs10000023          4  95733906     658 BMPR1B           + 0.8642906
rs10000030          4 103374154      NA   <NA>        <NA> 0.2832705 
}
}
\seealso{
\code{\link{genome_order}}
}
\examples{
## DEMO // WORKSPACE
data(demo,SNPsAnnotation)
all_data <- read_pvals(data_name=demo,snps_ann=SNPsAnnotation)

\dontrun{
##
##  Below is an example on how to annotate the SNPs prior the use of genomicper 
##  using UCSC table browser and intersectBed from bedtools:

##  The function intersectBed from bedtools can be used to annotate SNPs to genes. 
##  This function needs the locations to be annotated as input, and a reference file
##  to annotate to. Genomicper uses entrez gene ids to annotate associate SNPs-to genes-pathways. 
 
# prepare locations INPUT: chr	position	position	other-info

# 1       10763241        10763241        1_10763241_C_T_1
# 1       10764465        10764465        1_10764465_T_C_1
# 1       10767685        10767685        1_10767685_C_T_1

# Prepare the file to annotate to. Using UCSC table browser. 
# clade:Mammal	genome:Human	assembly: Feb2009(GRCh37/hg19)
# group: All tables	database:hg19 Table: knownToLocusLink 
# output format: selected fields from primary and related tables
# click on "get output" 
# Next select Linked Tables: kgXref and knownGene
# click on "allow filtering using fields in checked tables" 
# Select fields for output: 
	# Entrez Gene ID from hg19.knownToLocusLink 
	# Gene Symbol from hg19.kgXref
	# Reference sequence chromosome or scaffold from hg19.knownGene
	# + or - for strand from hg19.knownGene
	# Transcription start position from hg19.knownGene
	# Transcription end position from hg19.knownGene
# click on "get output"
# Table will include more than one mapping, to avoid results bias  decrease/increase 
# the min and max according to the wished annotations for a single gene 
# (eg. take min and max of all isoforms or desiered kb distance)

# Reformat Table to intersectBed accepted formats (eg.GTF/BED/VCF)
# awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$5,$6,$1,$2,$4}' Genes_hg19_TableBrowser.txt | 
# sed 's/chr//g' | awk 'BEGIN{FS="\t";OFS="\t"}{if($1 !~ /[:alnum:]/) print $0}' > Genes_TEMP.txt

# R > 
# x <- read.table("Genes_TEMP.txt",sep="\t",header=F,stringsAsFactors=F)
# genes <- unique(sort(x[,5]))
# gene_table <- matrix(data=NA,ncol=6,nrow=0)
# for(i in genes){
	# grids <- which(x[,5] == i)
	# min <- x[grids[which.min(x[grids,2])],2]
	# max <- x[grids[which.max(x[grids,3])],3]
	# gene_table <- rbind(gene_table,c(x[grids[1],1],min,max,
	#           x[grids[1],4],x[grids[1],5],x[grids[1],6])) 	
# }
# write.table(gene_table,file="Gene_Table.txt",col.names=F,row.names=F,sep="\t",quote=F)
# /exit R

## If you are trying to intersect very large files and are having trouble 
## with excessive memory usage, please presort your data by chromosome
## and then by start position e.g.: sort -k1,1 -k2,2n in.bed > in.sorted.bed 
## for BED files) and then use the -sorted option
## sort -k1,1 -k2,2n Gene_Table.txt > Gene_Table_sorted.txt

## Intersect command:
# intersectBed -a inp.txt -b Gene_Table_sorted.txt -wa -wb -sorted > temp
# Select Columns : SNP_ID,CHR,SNP_Location,GeneID,OtherAnnotation1,OtherAnnotation2
# awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$5,$2,$8,$9,$10}' temp > SNP_Table_Annotation.txt

# data ready for genomicper:
# head SNP_Table_Annotation.txt
# rs1000313       1       15405489        23254   KAZN    +
# rs1002365       1       19797248        832     CAPZB   -
# rs1002487       1       26865971        6195    RPS6KA1 +
# rs1002358       1       53753718        7804    LRP8    -
# rs1001160       1       76358591        4438    MSH4    +
# rs1002784       1       76824595        256435  ST6GALNAC3      +
# rs1001193       1       147166377       400818  NBPF9   +
# rs1001193       1       147166377       728841  NBPF8   +
# rs1001193       1       147166377       728855  LINC00623       +
# rs1001193       1       147166377       653505  PPIAL4B +
}
}
\keyword{ gwas }
\keyword{ pvalues }
