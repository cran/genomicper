get_results <-
function(res_pattern="Permus"){
permus_list <-NULL
y<-"NULL"
for(i in list.files(pattern="Permus")){
mode<-strsplit(i,split="[/.]")
mode<-mode[[1]]
##current trait
x <-strsplit(mode,split="[_]")[[1]]
x <-x[2]
x <- paste("Permus_",x,sep="")
mode<-mode[1]
#current trait the same?
if(x %in% permus_list){
z <- which(permus_list==x)
permus_list[z]
assign(permus_list[z],cbind(get(permus_list[z]),read.table(i)))
}
else{
assign(mode,read.table(i))
permus_list <- c(permus_list,mode)
}
}
ntraits <- length(permus_list)
temp <- get(permus_list[1])
npaths <- dim(temp)[2]
npermus <- dim(temp)[1]
top<-npermus+1
threshold <- 0.05
id_thre<-round(npermus*threshold)+1

reslts <- NULL
for(i in 1:dim(get(permus_list[1]))[2]){
reslts<- c(reslts,strsplit(colnames(get(permus_list[1])),"_")[[i]][3])
}
trts<-NULL
for(i in 1:length(permus_list)){
trts<-c(trts,strsplit(permus_list,"_")[[i]][2])
}

hyper_vs_empirical <- matrix(data=0,nrow=0,ncol=5)
colnames(hyper_vs_empirical) <- c("PathID","Trait","Threshold","P-Value","Observed")
for(i in 1:length(permus_list)){
temp <- get(permus_list[i])
mat <- matrix(data=0,nrow=npaths,ncol=5)
mat[,1] <- reslts
for(j in 1:npaths){
x <-  sort(temp[,j])
thre <- x[id_thre]
mat[j,2] <- trts[i]
mat[j,3] <- thre
pval <- temp[1,j]
mat[j,4] <- pval
indx <- which(x==pval)[1]
obs <- indx/npermus
mat[j,5] <- obs
}
hyper_vs_empirical<-rbind(hyper_vs_empirical,mat)
}

hyper_vs_empirical<- as.data.frame(hyper_vs_empirical,stringsAsFactors=FALSE)
hyper_vs_empirical[,3]<- as.numeric(hyper_vs_empirical[,3])
hyper_vs_empirical[,4]<- as.numeric(hyper_vs_empirical[,4])
hyper_vs_empirical[,5]<- as.numeric(hyper_vs_empirical[,5])
return(hyper_vs_empirical)
}

