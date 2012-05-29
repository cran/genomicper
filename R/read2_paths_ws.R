read2_paths_ws <-
function(ordered_alldata=ordered_alldata,gs_locs=gs_locs,pathinfo=pathways){
ps <- ls(pattern="hsa",envir=.GlobalEnv)
pathways <- pathinfo

ttl_paths <- length(ps)
paths_ids<- list()
paths_names <- NULL

per_ors <- list()
indx_paths<-NULL
per_orids_list <- NULL
pathways2<- matrix(data=0,ncol = 5,nrow=ttl_paths)


for(i in 1:ttl_paths){
moe <- ps[i]
paths_names <- c(paths_names,moe)
per_orids_list<-c(per_orids_list,paste("per_orids_",moe,sep=""))
genes <- get(ps[i])
paths_ids<- c(paths_ids,list(moe=genes))

indx_paths <-which(pathways[,1] == moe)
ids<-NULL
ids2<-NULL

pathways2[i,1:3] <- pathways[indx_paths,1:3]
for(j in 1:length(genes)){
x <- which(as.numeric(gs_locs[,4])== as.numeric(genes[j]))
if(length(x)!=0){
ids <- c(ids,x)
}
y <- which(as.numeric(as.character(ordered_alldata[,4]))== as.numeric(genes[j]))
if(length(y)!=0){
ids2 <- c(ids2,y)
}    
if(j == length(genes)){
pathways2[i,4]<-length(ids)
pathways2[i,5]<-length(ids2)
} 
}
per_ors <-c(per_ors,list(ids))
}
names(paths_ids) <- paths_names
names(per_ors) <- per_orids_list
colnames(pathways2)<-c("ID","Name","GenesInPath","GenesFound","SNPsInPath")
return(list(pathways2=pathways2,per_ors=per_ors))
}#close function

