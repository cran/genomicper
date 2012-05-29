genome_permutation <-
function(ordered_alldata="",pers_ids=paths_res$per_ors,pathways2=paths_res$pathways2,ntraits="",nper=100,threshold=0.05){
ntraits <- as.numeric(ntraits)
nper <- as.numeric(nper)
threshold <- as.numeric(threshold)
temp <- ordered_alldata[,c(1:6,ntraits)]
paths_list <- names(pers_ids)


## Set number and type of permutations:
# mx_rs holds the max number of rows in all_data ordered_alldata
mx_rs <- dim(ordered_alldata)[1]
set.seed(10, kind="Mersenne-Twister")
sd <-round(runif(nper,1,mx_rs))
rowsf <- dim(gs_locs)[1]

tname <- NULL
lab <- NULL
i <- NULL
ids<-NULL
j<-NULL
k<-NULL
date()
all_ts<-NULL
listf <- as.numeric(as.character(gs_locs[,4]))
pathways3 <- matrix(data=0,nrow=length(paths_list),ncol=5+length(ntraits))
pathways3[,1:5] <- pathways2[,1:5]
pathways2<-pathways3
rm(pathways3)

for(i in 7:length(temp)){ # open traits
hyper_mat<-matrix(data=NA,nrow=length(sd)+1,ncol=length(paths_list))
trtsnm <- colnames(temp)[i]
all_ts <- c(all_ts,trtsnm) 
sig_genes_real <- 0
for(j in 1:length(sd)){ # open permutations
sig_genes <- 0
path_counts <- matrix(data=0,nrow=length(paths_list),ncol=2)
for(m in 1:length(paths_list)){
path_counts[m,1]<- paths_list[m]
}
for(k in 1:length(listf)){ # open genes
obs <- as.numeric(gs_locs[k,6])
strt<- as.numeric(gs_locs[k,5])
ed <- strt + obs - 1
if(j == 1){
tscore <- pchisq(-2*(sum(log(temp[c(strt:ed),i]))),df=2*(obs),lower.tail=FALSE)
if(tscore <= threshold){
sig_genes_real <- sig_genes_real + 1
for(m in 1:length(paths_list)){
if(k %in% pers_ids[[m]]){
pathways2[m,i-1]<- as.numeric(pathways2[m,i-1])+1
}
}
}
}
obs <- as.numeric(gs_locs[k,6])
fkstrt <- as.numeric(gs_locs[k,5]) + sd[j]
fkend <- fkstrt + obs - 1
if(fkend <= mx_rs){ 
tscore <- pchisq(-2*(sum(log(temp[c(fkstrt:fkend),i]))),df=2*(obs),lower.tail=FALSE)
if(tscore <= threshold){
sig_genes <- sig_genes + 1
for(m in 1:length(paths_list)){
if(k %in% pers_ids[[m]]){
path_counts[m,2]<- as.numeric(path_counts[m,2])+1
}
}
next
}
}
if(fkend > mx_rs){ 
if(fkstrt >= mx_rs){ 
pls <- fkstrt - mx_rs 
nend <- pls + obs - 1
if(nend > mx_rs){ 
rst <- nend - mx_rs
ps  <- c(pls:mx_rs,1:rst)
tscore <- pchisq(-2*(sum(log(temp[ps,i]))),df=2*(obs),lower.tail=FALSE)
if(tscore <= threshold){
sig_genes <- sig_genes + 1
for(m in 1:length(paths_list)){
if(k %in% pers_ids[[m]]){
path_counts[m,2]<- as.numeric(path_counts[m,2])+1
}
}
next
}
} 
else{ 
ps <- c(pls:nend)
tscore <- pchisq(-2*(sum(log(temp[ps,i]))),df=2*(obs),lower.tail=FALSE)
if(tscore <= threshold){
sig_genes <- sig_genes + 1
for(m in 1:length(paths_list)){
if(k %in% pers_ids[[m]]){
path_counts[m,2]<- as.numeric(path_counts[m,2])+1
}
}
next
}
} 
} 
else{ 
rst <- fkend - mx_rs
ps  <- c(fkstrt:mx_rs,1:rst)
tscore <- pchisq(-2*(sum(log(temp[ps,i]))),df=2*(obs),lower.tail=FALSE)
if(tscore <= threshold){
sig_genes <- sig_genes + 1
for(m in 1:length(paths_list)){
if(k %in% pers_ids[[m]]){
path_counts[m,2]<- as.numeric(path_counts[m,2])+1
}
}
}
} 
} 
}#close genes
if(j==1){
for(m in 1:length(paths_list)){
hypr_ps <- HyperBigUniverse_v2(Sig_in_Paths=as.numeric(pathways2[m,i-1]),uniSig=sig_genes_real,gns_in_Paths=as.numeric(pathways2[m,4]),universe=rowsf)
hyper_mat[j,m]<- hypr_ps
hypr_ps <- HyperBigUniverse_v2(Sig_in_Paths=as.numeric(path_counts[m,2]),uniSig=sig_genes,gns_in_Paths=as.numeric(pathways2[m,4]),universe=rowsf)
hyper_mat[j+1,m]<- hypr_ps
}
}
else{
for(m in 1:length(paths_list)){
hypr_ps <- HyperBigUniverse_v2(Sig_in_Paths=as.numeric(path_counts[m,2]),uniSig=sig_genes,gns_in_Paths=as.numeric(pathways2[m,4]),universe=rowsf)
hyper_mat[j+1,m]<-hypr_ps
}
}
}# close premutations
rownames(hyper_mat)<-c("REAL_Pval",1:length(sd))
colnames(hyper_mat)<- paths_list
write.table(hyper_mat,file=paste("Permus_",trtsnm,".txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
}#close traits
colnames(pathways2)<-c("ID","Name","GenesInPath","GenesFound","SNPsInPath",all_ts)
return()
}

