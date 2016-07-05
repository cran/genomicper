get_pathways <-
function(source="reactome",all_paths=TRUE,envir = ""){
if(missing(source)==FALSE){
if(source != "reactome" && source != "kegg"){
stop("Argument source must be set to:\"reactome\" OR \"kegg\"")
}
}
print("Arguments set:")
    print(paste("Source:" ,substitute(source),sep=""))
    print(paste("Annotate all pathways: ",substitute(all_paths),sep=""))
if(source=="reactome"){
pathways_description <- dbGetQuery(reactome.db::reactome_dbconn(), "SELECT * FROM pathway2name")
x<-unique(sub(":.*","",pathways_description[,2]))
print(x) #"select species"
sp <- readline(prompt="Select dataset(e.g. \"Homo sapiens\"):\n")
sp <- gsub("\"","",sp)
x<- grep(sp, pathways_description[,2],ignore.case=TRUE)
pathways_description <- pathways_description[x,] 
pathways_description[,2]<- sub(".*: ","",pathways_description[,2])
paths <- AnnotationDbi::as.list(reactome.db::reactomePATHID2EXTID)
pathways_description <- pathways_description[which(pathways_description[,1]%in% names(paths) == TRUE),]
pre <- readline(prompt="Assign a Prefix for the pathways (e.g.\"mypath\",\"rhsa\"):\n")
if(all_paths==TRUE){
tt <- dim(pathways_description)[1]
print(paste(tt," reactome pathways to annotate",sep=""))
for(i in 1:tt){
print(pathways_description[i,2])
x <- which(names(paths)== pathways_description[i,1])
if(length(x)!=0){
assign(paste(pre,pathways_description[i,1],sep=""),paths[[x]], envir = envir)
}
}
}
if(all_paths==FALSE){
print(pathways_description)
subs <- readline(prompt="Select reactome IDs - comma separated (e.g. \"70153,69618\"):\n")
subs <- strsplit(subs,",")[[1]]
subs <- gsub("\"","",subs)
tt <- length(subs)
y <- NULL
z <- NULL
for(i in 1:tt){
exs <- which(pathways_description[,1]== subs[i])
if(length(exs)!=0){
print(pathways_description[exs,2])
x <- which(names(paths)== pathways_description[exs,1])
if(length(x)!=0){
assign(paste(pre,pathways_description[exs,1],sep=""),paths[[x]], envir = envir)
z <- c(z,exs)
}

}
else{
y <- c(subs[i],y)
}
}
print(paste("Pathways without gene ids/not found:",y,sep=""))
pathways_description <- pathways_description[z,]
}
print(paste("Pathwas saved in workspace as: ",pre,"XXXXXX",sep=""))
}
if(source=="kegg"){
paths <- AnnotationDbi::as.list(KEGG.db::KEGGPATHID2EXTID)
x <- names(paths)
print(unique(substr(x,1,3))) #"select species"
sp <- readline(prompt="Select the organism code: (e.g. \"hsa\"):\n")
sp <- gsub("\"","",sp)
ids <- grep(sp,x,ignore.case = TRUE)
paths <- paths[ids] 
ps<-length(paths)
path_names<-AnnotationDbi::as.list(KEGG.db::KEGGPATHID2NAME)
pathways_description <- matrix(data=NA, nrow=ps, ncol=2)
for(i in 1:ps){
pathways_description[i,1]<- names(paths)[i]
pathways_description[i,2]<- path_names[[which(names(path_names)== substr(names(paths)[i],4,8))]]
}
if(all_paths==TRUE){
tt <- length(ids)
print(paste(tt," reactome pathways to annotate",sep=""))
for(i in 1:tt){
print(names(paths)[i])
assign(names(paths)[i],paths[[i]], envir = envir)
}
}
if(all_paths==FALSE){
print(pathways_description)
subs <- readline(prompt="Select pathways comma separated (e.g. \"hsa00232,hsa00230\"):\n")
subs <- strsplit(subs,",")[[1]]
subs <- gsub("\"","",subs)
tt <- length(subs)
y <- NULL
z <- NULL
for(i in 1:tt){
x <- which(names(paths)== subs[i])
if(length(x)!=0){
print(names(paths)[x])
assign(names(paths)[x],paths[[x]], envir = envir)
z <- c(z,x)
}
else{
y <- c(subs[i],y)
}
}
print(paste("Pathways not found:",y,sep=""))
pathways_description <- pathways_description[z,]
}
print(paste("Pathwas saved in workspace as: ",sp,"XXXXXX",sep=""))
}
return(pathways_description)
}
