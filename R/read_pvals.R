read_pvals <-
function(data_name="",snps_ann="",snps_ann_all=""){
data <- read.table(data_name,sep="\t",header=T)
all_snps <- read.table(snps_ann,sep="\t",header=T)
all_snps_update <- read.table(snps_ann_all,sep="\t",header=T)

# Set the names of the first column (SNP_IDS) to be the same 'name' 
colnames(data)[1]<- "name"
colnames(all_snps)[1]<-"name"
colnames(all_snps_update)[1] <-"name"
# Merge GWAS p-values and SNPs to gene Annotations at Distance 0
# Note: Number of SNPs may increase, due to SNPs annotated to more than one gene)
all_data<-merge(all_snps,data, by="name",all.x=FALSE,all.y=TRUE, incomparables="NA")

# Merge SNPs to gene Annotations at Distance 0 with ALL SNPS annotations to update the non-annotated-SNP to genes
all_data2<-merge(all_data,all_snps_update, by="name",all.x=TRUE,all.y=FALSE, incomparables="NA")

## Update annotations 
#last two columns: chromosome and Location from the "SNPs2Genes_ALL.txt"
chrms <- dim(all_data2)[2]-1
locs <- dim(all_data2)[2]
# Replace with the new ones
all_data2[,2]<-all_data2[,chrms]
all_data2[,3]<-all_data2[,locs]
# Subset data 
cols <- dim(all_data)[2]
all_data2 <- all_data2[,1:cols]

# Number of SNPS merged or modified to other SNPs
print("Number of SNPS without genomic location:")
print(length(which(is.na(all_data2[,2])== TRUE )))
all_data <- all_data2
return(all_data)
}

