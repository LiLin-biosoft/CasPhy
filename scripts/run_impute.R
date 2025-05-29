suppressMessages(library(CasPhy))
suppressMessages(library(parallel))
args<-commandArgs(trailingOnly=T)

fun1<-function(i){
  findNeighborUnweightedParallel(
    i=i-1,mat=binary,maxChar=maxChar,
    n=n_neighbor,p=p_mutation
  )
}
fun2<-function(i){
  findNeighborRangeUnweighted(
    i=i-1,mat=binary,maxChar=maxChar,
    range=range_tmp,offset=offset,
    n=n_neighbor,p=p_mutation
  )
}
impute.fun<-function(i){
  imputeMat(i-1,mat,neighbor_found,offset)
}

mat<-as.matrix(read.table(file=args[1]))
rownames(mat)<-paste0("C",1:nrow(mat))
binary<-mat2bit(mat)$binary
output<-matrix(NA,nrow=nrow(mat),ncol=ncol(mat))
maxChar<-ncol(mat)

status<-read.table(file=args[3],row.names=1)
run_mode<-status["run_mode",1]
n_iter<-status["n_iter",1]
n_neighbor<-status["n_neighbor",1]
p_mutation<-status["p_mutation",1]
batch_size<-status["batch_size",1]
ncore<-as.integer(args[4])

cat(as.character(Sys.time()),": iteration ",n_iter," start, n_neighbor = ",n_neighbor,sep="",fill=T)
if(run_mode==2){
  range.desc<-NULL
  for(i in 1:ceiling(nrow(binary)/batch_size)-1){
    cat(as.character(Sys.time()),": batch ",i+1," start",sep="",fill=T)
    offset<-i*batch_size
    row_process<-(offset+1):min(offset+batch_size,nrow(binary))
    cl<-makeForkCluster(ncore)
    neighbor_found<-t(parSapply(cl,row_process,fun1))
    stopCluster(cl)
    cat(as.character(Sys.time()),": finding neighbor finished",sep="",fill=T)
    cl<-makeForkCluster(ncore)
    output[row_process,]<-t(parSapply(cl,row_process,impute.fun))
    stopCluster(cl)
    cat(as.character(Sys.time()),": imputing mat finished",sep="",fill=T)
    rm(neighbor_found)
  }
  temp<-sum(mat==-1)-sum(output==-1)
  cat("Impute ",temp," entries.",sep="",fill=T)
  write.table(output,file=args[2],row.names=F,col.names=F,quote=F,sep="\t")
  status["run_mode",1]<-1
  status["n_iter",1]<-n_iter+1
  if(temp<max(floor(nrow(mat)/100),100)){
    status["n_neighbor",1]<-n_neighbor-1
  }
  write.table(status,file=args[3],row.names=T,col.names=F,quote=F,sep="\t")
}else if(run_mode==1){
  full.bk<-paste0(args[5],".bk")
  full.desc<-paste0(args[5],".desc")
  if(file.exists(full.bk)){
    unlink(full.bk)
  }
  if(file.exists(full.desc)){
    unlink(full.desc)
  }
  filepath<-gsub("[^/]+$","",args[5])
  filename<-gsub(filepath,"",args[5])
  range.bk<-big.matrix(
    nrow=ceiling(nrow(mat)/32),
    ncol=nrow(mat),
    type="integer",
    backingfile=paste0(filename,".bk"),
    descriptorfile=paste0(filename,".desc"),
    backingpath=filepath
  )
  range.desc<-describe(range.bk)
  cat("initialing",fill=T)
  for(i in 1:ceiling(nrow(binary)/batch_size)-1){
    cat(as.character(Sys.time()),": batch ",i+1," start",sep="",fill=T)
    offset<-i*batch_size
    row_process<-(offset+1):min(offset+batch_size,nrow(binary))
    cl<-makeForkCluster(ncore)
    neighbor_found<-t(parSapply(cl,row_process,fun1))
    stopCluster(cl)
    cat(as.character(Sys.time()),": finding neighbor finished",sep="",fill=T)
    cl<-makeForkCluster(ncore)
    output[row_process,]<-t(parSapply(cl,row_process,impute.fun))
    stopCluster(cl)
    cat(as.character(Sys.time()),": imputing mat finished",sep="",fill=T)
    range.bk[,row_process]<-t(neighbor_found)
    cat(as.character(Sys.time()),": writing range finished",sep="",fill=T)
    rm(neighbor_found)
  }
  temp<-sum(mat==-1)-sum(output==-1)
  cat("Impute ",temp," entries.",fill=T)
  write.table(output,file=args[2],row.names=F,col.names=F,quote=F,sep="\t")
  status["run_mode",1]<-0
  status["n_iter",1]<-n_iter+1
  if(temp<max(floor(nrow(mat)/100),100)){
    status["n_neighbor",1]<-n_neighbor-1
  }
  write.table(status,file=args[3],row.names=T,col.names=F,quote=F,sep="\t")
}else{
  range.bk<-attach.big.matrix(
    paste0(args[5],".desc")
  )
  range.desc<-describe(range.bk)
  for(i in 1:ceiling(nrow(binary)/batch_size)-1){
    cat(as.character(Sys.time()),": batch ",i+1," start",sep="",fill=T)
    offset<-i*batch_size
    row_process<-(offset+1):min(offset+batch_size,nrow(binary))
    range_tmp<-t(range.bk[,row_process])
    cat(as.character(Sys.time()),": reading range finished",sep="",fill=T)
    cl<-makeForkCluster(ncore)
    neighbor_found<-t(parSapply(cl,row_process,fun2))
    stopCluster(cl)
    cat(as.character(Sys.time()),": finding neighbor finished",sep="",fill=T)
    cl<-makeForkCluster(ncore)
    output[row_process,]<-t(parSapply(cl,row_process,impute.fun))
    stopCluster(cl)
    cat(as.character(Sys.time()),": imputing mat finished",sep="",fill=T)
    rm(neighbor_found,range_tmp)
  }
  temp<-sum(mat==-1)-sum(output==-1)
  cat("Impute ",temp," entries.",sep="",fill=T)
  write.table(output,file=args[2],row.names=F,col.names=F,quote=F,sep="\t")
  status["n_iter",1]<-n_iter+1
  if(temp<max(floor(nrow(mat)/100),100)){
    if(n_neighbor==5){
      status["run_mode",1]<-1
    }
    status["n_neighbor",1]<-n_neighbor-1
  }
  write.table(status,file=args[3],row.names=T,col.names=F,quote=F,sep="\t")
}

