get_results_ws <-
function(res_pattern="Permus"){
permus_list <-NULL
y<-"NULL"
permus_list <- ls(pattern="Permus",envir=.GlobalEnv)
ntraits <- length(permus_list)
print(permus_list)
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

