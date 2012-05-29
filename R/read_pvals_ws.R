read_pvals_ws <-
function(data_name="",snps_ann="",snps_ann_all=""){
data <- data_name
all_snps <- snps_ann
all_snps_update <- snps_ann_all
all_data<-merge(all_snps,data, by="name",all.x=FALSE,all.y=TRUE, incomparables="NA")
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

