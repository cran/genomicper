genome_order<-
function(all_data=""){
# Number of unique genes
rowsf<-length(unique(sort(as.character(all_data[,4]))))
# Total number of columns 
colsf<-length(all_data)
## Prepare data
or_data <- all_data
or_data[,2] <- as.character(or_data[,2])
xs<-which(or_data[,2]=="X")
if(length(xs)!=0){
or_data[xs,2]<-"23"
}
ys<-which(or_data[,2]=="Y")
if(length(ys)!=0){
or_data[ys,2]<-"24"
}
or_data[,2] <- as.numeric(or_data[,2])

nas <- or_data[,3]
nas[is.na(nas)] <- 0
or_data[,3]<-as.numeric(nas)

print(paste("Number of SNPs Removed:",length(which(or_data[,3]== "0")),sep=""))

if(class(or_data[,2])!="numeric"){
or_data[,2]<- as.numeric(as.character(or_data[,2]))
}

if(class(or_data[,3])!="numeric"){
or_data[,3]<- as.numeric(as.character(or_data[,3]))
}

## Order data (1) to take the smallest number on the location of the genes
indx <- order(or_data[,2],or_data[,3])
or_data<- or_data[indx,]

## unique genes constant
listf <- unique(sort(or_data[,4]))
### Get Gene Start Indexes and number of observations
# matrix of obs
gs_locs <- matrix(data=NA,nrow=rowsf,ncol=4)

for(j in 1:length(listf)){
ids <- which(or_data[,4]==listf[j])
obs <- length(ids)
gs_locs[j,1]<- as.character(or_data[ids[1],5])
gs_locs[j,2]<- or_data[ids[1],2]
gs_locs[j,3]<- or_data[ids[1],3]
gs_locs[j,4]<- listf[j]
}

lab<-c("Symbol","Chromosome","Location1st","GENE_ID")
colnames(gs_locs)<-lab

colnames(gs_locs)
ordered_alldata<-merge(gs_locs,or_data, by="GENE_ID",all=TRUE, incomparables="NA")

## Change to numeric
ordered_alldata[,3]<-as.numeric(as.character(ordered_alldata[,3]))
ordered_alldata[,6]<-as.numeric(as.character(ordered_alldata[,6]))
ordered_alldata[,4]<-as.numeric(as.character(ordered_alldata[,4]))
ordered_alldata[,7]<-as.numeric(as.character(ordered_alldata[,7]))

### COPY SNP LOCATION TO THOSE WITHOUT GENE ANNOTATIONS
strt<- which(is.na(ordered_alldata[,4])==TRUE)[1]
endf<-dim(ordered_alldata)[1]

## copy
ordered_alldata[strt:endf,4]<- ordered_alldata[strt:endf,7]

## Re-ordering data has to be using the Locations1st column; must be numeric (column 4)
## First Order by Chromosome then by Location_1st
#### sorting 1st the chromosome, 2nd the location
indx <- order(ordered_alldata[,6],ordered_alldata[,4])
ordered_alldata2<- ordered_alldata[indx,]

# keep (name,Chromosome.y,Locations1st, GENE_ID,Symbol.x,Orientation,...traits)
ordered_alldata2<- ordered_alldata2[,c(5,6,4,1,2,9,10:length(ordered_alldata))]
names(ordered_alldata2)

colnames(ordered_alldata2)[2]<-"Chromosome"
colnames(ordered_alldata2)[5]<-"Symbol"

# Remove SNPs merged wth the update/ SNPs without genomic locations
mx_rs <- which(is.na(ordered_alldata2[,2])== TRUE)[1]-1
if(is.na(mx_rs)==FALSE){
ordered_alldata2 <- ordered_alldata2[1:mx_rs,]
rownames(ordered_alldata2)<- c(1:mx_rs)
}
if(is.na(mx_rs)==TRUE){
endf2 <- dim(ordered_alldata2)[1]
rownames(ordered_alldata2)<- c(1:endf2)
}
###################################
### Make counts second time to find starting indexes for each gene and # of observations
gs_locs <- matrix(data=NA,nrow=rowsf,ncol=6)

for(j in 1:length(listf)){
ids <- which(ordered_alldata2[,4]== listf[j])
obs <- length(ids)
gs_locs[j,1]<- as.character(ordered_alldata2[ids[1],5])
gs_locs[j,2]<- ordered_alldata2[ids[1],2]
gs_locs[j,3]<- ordered_alldata2[ids[1],3]
gs_locs[j,4]<- listf[j]
gs_locs[j,5]<- ids[1]
gs_locs[j,6]<- obs
}

lab<-c("Symbol","Chromosome","Location","Gene_ID","Start_Indx","Observations")
colnames(gs_locs)<-lab
## gs_locs[,5] = START ID from ordered_alldata
## gs_locs[,6] = Number of observatoins per gene for ordered_alldata
########################################################################
ordered_alldata<-ordered_alldata2
return(list(ordered_alldata=ordered_alldata,gs_locs=gs_locs))
}

