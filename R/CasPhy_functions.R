colsfun<-function(n,seed=NA){
  pal<-c("red","chocolate","darkorange","gold","darkgreen","green3","darkolivegreen1",
         "blue","lightblue","cyan","lightpink","magenta","purple")
  cols<-colorRampPalette(pal)(n)
  if(is.na(seed))
    return(cols)
  else{
    set.seed(seed)
    return(cols[sample(1:length(cols))])
  }
}

gen.rowCol<-function(splt,cols){
  temp<-as.integer(factor(splt[,1],levels=sort(unique(c(0,splt[,1])))))
  temp[temp>length(cols)]<-length(cols)+1
  output<-c("white",cols,"grey")[temp]
  for(i in 2:ncol(splt)){
    output<-cbind(output,"white","white")
    for(j in sort(setdiff(unique(splt[,i-1]),0))){
      idx<-which(splt[,i-1]==j)
      temp<-splt[idx,i]
      temp<-as.integer(factor(temp),levels=sort(unique(c(0,temp))))
      temp[temp>length(cols)]<-length(cols)+1
      output[idx,i*2-1]<-c("white",cols,"grey")[temp]
    }
  }
  return(output)
}

splt2col<-function(splt,cols){
  if(is.vector(splt)){
    splt<-matrix(splt,,1)
  }
  odr<-do.call(order,as.data.frame(splt))
  clust_cols<-splt
  clust_cols[clust_cols!=0]<-clust_cols[clust_cols!=0]%%2+1
  clust_cols<-clust_cols+1
  clust_cols<-matrix(cols[clust_cols],nrow(splt),ncol(splt))
  return(list(odr=odr,cols=clust_cols))
}

#' Turns the matrix into bit-coded data.
#'
#' @param mat The mutation matrix.
#' @param code_unmutated Whether to code unmutated state 0.
#' @param return.bit Whether to return the bit-coded matrix.
#' 
#' @examples
#' \dontrun{
#' binary<-mat2bit(mat)
#'}
#'
#' @return A list with the following components:
#' @return - Binary: is the bit-coded data.
#' @return - loci: is the offset of each locus.
#' @return - char: is the symbols in the original mutation matrix.

mat2bit<-function(mat,code_unmutated=F,return.bit=T){
  if(code_unmutated){
    fixed_char<-c(-1)
  }else{
    fixed_char<-c(-1,0)
  }
  char<-lapply(1:ncol(mat),function(i){
    sort(setdiff(unique(mat[,i]),fixed_char))
  })
  loci<-diffinv(sapply(char,length))
  if(return.bit){
    for(i in 1:ncol(mat)){
      tmp<-c(fixed_char,char[[i]])
      if(code_unmutated){
        mat[,i]<-match(mat[,i],tmp)-1
      }else{
        mat[,i]<-match(mat[,i],tmp)-2
      }
    }
    binary<-matCoding(mat,loci)
    rownames(binary)<-rownames(mat)
    return(list(binary=binary,loci=loci,char=unname(unlist(char))))
  }else{
    return(list(binary=NULL,loci=loci,char=unname(unlist(char))))
  }
}

regularize.split<-function(x,min_size=0){
  category<-unique(x)
  x<-factor(x,levels=category)
  levels(x)<-match(1:nlevels(x),order(table(x),decreasing=T))
  x<-as.integer(as.character(x))
  if(min_size>1){
    idx<-(1:max(x))[table(x)<min_size]
    x[x%in%idx]<-0
  }
  return(x)
}

#' The top solver that splits the data into small clusters
#'
#' @param mat The mutation matrix.
#' @param bit_coded Whether the matrix is bit coded.
#' @param weight Weights of sites.
#' @param splt The results of `split.data` from previous run.
#' @param split_idx The index of cluster to further split.
#' @param radius The minimum number of mutations between neighbors.
#' @param p_shared The threshold for the percentage of shared neighbors.
#' @param min_size The minimum size of clusters. Clusters smaller than this value will be collapsed.
#' @param file_backed The name of backup file. If set to `NULL`, do not create backup file. If set to `""`, create an anonymous backup file.
#' @param ncore The number of threads used.
#' @param cuda_path The path of CUDA script. If set to `NULL`, do not use CUDA.
#' @param verbose Whether print detailed information.
#' 
#' @examples
#' \dontrun{
#' splt<-split.data(mat,radius=25,p_shared=0.8)
#'}
#'
#' @return Return an integer matrix containing the cluster information.

split.data<-function(mat,bit_coded=F,weight=NULL,splt=NULL,splt_idx=1,
                     radius=25,p_shared=0.8,min_size=10,
                     file_backed=NULL,ncore=1,cuda_path=NULL,verbose=F){
  batch_size<-1e4
  if(!bit_coded){
    binary<-mat2bit(mat)
    if(!is.null(weight)&class(weight)[1]=="matrix"){
      temp<-diff(binary$loci)
      weight<-unlist(lapply(which(temp!=0),function(i){
        weight[1:temp[i],i]
      }))
    }
    binary<-binary$binary
    maxChar<-ncol(mat)
  }else{
    binary<-mat
    maxChar<-ncol(binary)*32
  }
  rm(mat)
  if(is.null(weight)){
    weight<-rep(1,ncol(binary)*32)
  }
  
  if(is.null(splt)){
    layer<-1
    k<-1
    splt<-cbind(rep(1,nrow(binary)),rep(0,nrow(binary)))
  }else{
    layer<-ncol(splt)
    layer<-which(apply(splt,2,function(x){
      splt_idx%in%x
    }))
    k<-max(splt)
    if(layer==ncol(splt)){
      splt<-cbind(splt,rep(0,nrow(binary)))
    }
  }
  
  idx<-which(splt[,layer]==splt_idx)
  binary<-binary[idx,]
  
  if(is.null(cuda_path)){
    neighbor.fun<-function(i){
      radiusWeightedParallel(
        i=i-1,mat=binary,maxChar=maxChar,
        weight=weight,radius=radius
      )
    }
    if(!is.null(file_backed)){
      if(file_backed==""){
        neighbor_mat<-big.matrix(
          nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),
          type="integer",backingfile=""
        )
        backfile<-file.name(neighbor_mat)
      }else{
        fpath<-gsub("[^/]+$","",file_backed)
        fname<-paste0(gsub(fpath,"",file_backed),".bk")
        neighbor_mat<-big.matrix(
          nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),type="integer",
          backingfile=fname,backingpath=fpath,descriptorfile=paste0(fname,".desc")
        )
        backfile<-paste0(file_backed,".bk")
      }
      if(verbose){
        cat("backup file:",backfile,fill=T)
      }
    }else{
      neighbor_mat<-big.matrix(nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),type="integer",backingfile=NULL)
    }
    nNeighbor<-integer()
    for(i in 1:ceiling(nrow(binary)/batch_size)-1){
      offset<-i*batch_size
      row_process<-(offset+1):min(offset+batch_size,nrow(binary))
      if(ncore>1){
        cl<-makeForkCluster(ncore)
        neighbor_tmp<-parSapply(cl,row_process,neighbor.fun)
        stopCluster(cl)
        rm(cl)
      }else{
        neighbor_tmp<-sapply(row_process,neighbor.fun)
      }
      neighbor_mat[,row_process]<-neighbor_tmp[-1,]
      nNeighbor<-c(nNeighbor,neighbor_tmp[1,])
      rm(neighbor_tmp)
      gc()
      if(verbose){
        cat(as.character(Sys.time()),": batch",i+1,"finished",fill=T)
      }
    }
    rm(binary,weight)
    gc()
    ncoreMerge<-min(ncore,ceiling(length(idx)/2e4))
    if(verbose){
      cat(as.character(Sys.time()),": start tipMerge with ",ncoreMerge," threads.",sep="",fill=T)
    }
    if(length(idx)<=batch_size||ncoreMerge==1||ncoreMerge==2){
      temp<-tipMergeNeighborBig(neighbor_mat@address,nNeighbor=nNeighbor,
                                p=p_shared,sqrt_nNeighbor=F,verbose=verbose)
    }else{
      nbatch<-ncoreMerge
      bsize<-ceiling(ncol(neighbor_mat)/nbatch)
      desc<-describe(neighbor_mat)
      cl<-makeForkCluster(ncoreMerge)
      mb<-parLapply(cl,1:nbatch-1,function(i){
        bgmat<-attach.big.matrix(desc)
        tipMergeSub(
          bgmat@address,
          i*bsize,bsize,
          nNeighbor,
          p=p_shared
        )
      })
      stopCluster(cl)
      if(verbose){
        cat("=")
      }
      mb<-member.combine(mb)
      tmp<-matrix(,nbatch,nbatch)
      colIdx<-col(tmp)[which(lower.tri(tmp))]-1
      rowIdx<-row(tmp)[which(lower.tri(tmp))]-1
      trans<-0:max(mb)
      offset<-0
      ntotal<-(nbatch*(nbatch-1)/2)
      while(offset<ntotal){
        processIdx<-(offset+1):min(offset+nbatch,ntotal)
        cl<-makeForkCluster(ncoreMerge)
        transTmp<-parSapply(cl,processIdx,function(i){
          bgmat<-attach.big.matrix(desc)
          offset1<-(colIdx[i])*bsize
          offset2<-(rowIdx[i])*bsize
          res<-tipMergeCross(
            bgmat@address,
            offset1,offset2,bsize,
            nNeighbor,mb,
            trans,p=p_shared
          )
          return(res)
        })
        stopCluster(cl)
        while(ncol(transTmp)>10){
          nRes<-ncol(transTmp)
          cl<-makeForkCluster(floor(nRes/2))
          if(nRes%%2==0){
            res<-parApply(cl,matrix(transTmp,,nRes/2),2,function(x){
              x<-matrix(x,,2)
              return(combineTransition(x[,1],x[,2]))
            })
            transTmp<-res
          }else{
            res<-parApply(cl,matrix(transTmp[,1:(nRes-1)],,(nRes-1)/2),2,function(x){
              x<-matrix(x,,2)
              return(combineTransition(x[,1],x[,2]))
            })
            transTmp<-cbind(res,transTmp[,nRes])
          }
          stopCluster(cl)
        }
        trans<-transTmp[,1]
        for(i in 2:ncol(transTmp)){
          trans<-combineTransition(trans,transTmp[,i])
        }
        offset<-offset+nbatch
        if(verbose){
          cat("=")        
        }
      }
      if(verbose){
        cat("",fill=T)
      }
      temp<-sapply(mb,function(i){
        transitionRoot(trans,i)
      })
    }
    rm(neighbor_mat)
    gc()
    if(verbose){
      cat(as.character(Sys.time()),": tipMerge finished.",sep="",fill=T)
    }
    if(!is.null(file_backed)){
      if(file.exists(backfile)){
        unlink(backfile)
      }
      if(file.exists(paste0(backfile,".desc"))){
        unlink(paste0(backfile,".desc"))
      }
    }
  }else{
    if(is.null(file_backed)||file_backed==""){
      cat("file_backed must be specified when using cuda.",fill=T)
      return(NULL)
    }
    write.mat2bin(binary,file=paste0(file_backed,".bmat"))
    writeBin(object=as.vector(weight),con=paste0(file_backed,".w"))
    rm(binary)
    cat(as.character(Sys.time()),": CUDA start.",sep="",fill=T)
    system(paste0(cuda_path,"radiusNeighbor ",file_backed," ",radius))
    system(paste0(cuda_path,"tipMerge ",file_backed," ",p_shared))
    cat(as.character(Sys.time()),": CUDA stop.",sep="",fill=T)
    unlink(paste0(file_backed,".bmat"))
    unlink(paste0(file_backed,".w"))
    unlink(paste0(file_backed,".bk"))
    unlink(paste0(file_backed,".nNeighbor"))
    unlink(paste0(file_backed,".idx"))
    temp<-read.table(file=paste0(file_backed,".member"))[,1]
    unlink(paste0(file_backed,".member"))
  }
  temp<-regularize.split(temp,min_size=min_size)
  temp[temp!=0]<-temp[temp!=0]+k
  splt[idx,layer+1]<-temp
  return(splt)
}

split.mg<-function(binary,maxChar,weight,n=50,p_mutation=0.5,radius=0,p_shared=0.8,sqrt_nNeighbor=F,min_size=10,verbose=F){
  if(is.null(weight)){
    neighbor_mat<-findNeighborUnweighted(
      binary,maxChar=maxChar,
      n=n,p=p_mutation,radius=radius,verbose=verbose
    )
  }else{
    neighbor_mat<-findNeighborWeighted(
      binary,maxChar=maxChar,weight=weight,
      n=n,p=p_mutation,radius=radius,verbose=verbose
    )
  }
  temp<-tipMergeNeighbor(neighbor_mat,p=p_shared,sqrt_nNeighbor=sqrt_nNeighbor,verbose=verbose)
  temp<-regularize.split(temp,min_size=min_size)
  return(temp)
}

#' impute matrix according to a phylogeny.
#'
#' @param mat The mutation matrix.
#' @param phy The phylogeny.
#' @param return.ans Whether return the ancestral states of nodes of the phylogeny.
#' 
#' @examples
#' \dontrun{
#' imputed<-phy.impute(mat,tre)
#'}
#'
#' @return an ultrametric phylogeny.

phy.impute<-function(mat,phy=NULL,return.ans=F){
  if(is.null(phy)){
    output<-apply(mat,2,function(x){
      temp<-x[x!=-1]
      if(length(unique(temp))==1){
        x<-rep(temp[1],length(x))
      }
      return(x)
    })
    rownames(output)<-rownames(mat)
    return(output)
  }
  mat_node<-node.ancestral(phy,mat)
  mat_total<-rbind(mat,mat_node)
  for(i in 1:phy$Nnode+length(phy$tip.label)){
    temp<-phy$edge[phy$edge[,1]==i,2]
    mat_total[temp,]<-t(apply(mat_total[temp,],1,function(x){
      imputeVector(ref=mat_total[i,],query=x)
    }))
  }
  if(return.ans){
    return(mat_total)
  }else{
    return(mat_total[1:nrow(mat),])
  }
}

#' turn the phylogeny into ultrametric phylogeny.
#'
#' @param phy The phylogeny object.
#' 
#' @examples
#' \dontrun{
#' tre<-as.ultrametric(tre)
#'}
#'
#' @return an ultrametric phylogeny.

as.ultrametric<-function(phy){
  e1<-phy$edge[,1]
  e2<-phy$edge[,2]
  if(is.null(phy$edge.length)){
    phy$edge.length<-rep(1,Nedge(phy))
  }else{
    phy$edge.length[phy$edge.length==0]<-phy$edge.length[phy$edge.length==0]+0
  }
  EL<-phy$edge.length
  depth<-numeric(length(phy$tip.label)+phy$Nnode)
  for (i in seq_len(length(e1))){
    depth[e2[i]]<-depth[e1[i]]+EL[i]
  }
  max_depth<-max(depth)
  delta_len<-apply(phy$edge,1,function(x){
    if(x[2]<=length(phy$tip.label)){
      return(max_depth-depth[x[2]])
    }else{
      return(0)
    }
  })
  phy$edge.length<-phy$edge.length+delta_len
  return(phy)
}

node.ancestral<-function(phy,mat){
  output<-matrix(,length(phy$tip.label)+phy$Nnode,ncol(mat))
  output[1:length(phy$tip.label),]<-mat
  for(i in phy$Nnode:1+length(phy$tip.label)){
    idx<-phy$edge[phy$edge[,1]==i,2]
    output[i,]<-charAncestral(output[idx,])
  }
  return(output[1:phy$Nnode+length(phy$tip.label),])
}

#' relabel the indels.
#'
#' @param mat The input matrix.
#' @param ref the reference matrix containg all indels.
#' 
#' @examples
#' \dontrun{
#' mat<-regularize.mat(mat,mat)
#'}
#'
#' @return an integer matrix.

regularize.mat<-function(mat,ref){
  for(i in 1:ncol(mat)){
    tmp<-sort(unique(c(ref[,i],-1L,0L)))
    mat[,i]<-match(mat[,i],tmp)-2L
  }
  return(mat)
}

# regularize.mat<-function(mat,ref){
#   for(i in 1:ncol(mat)){
#     tmp<-sort(unique(c(ref[,i],-1,0)))
#     mat[,i]<-match(mat[,i],tmp)-2
#   }
#   return(mat)
# }

greedy.wrap<-function(mat,weight=NULL,noise=F,collapse=F,file_backed=NULL,verbose=F){
  if(is.null(weight)){
    weight<-matrix(1,max(mat),ncol(mat))
  }
  if(noise){
    weight<-weight+matrix(runif(nrow(weight)*ncol(weight))/1e4,nrow(weight),ncol(weight))
  }
  phy<-list()
  if(!is.null(file_backed)){
    if(file_backed==""){
      shared_mat<-big.matrix(
        nrow=nrow(mat)*2-1,ncol=nrow(mat)*2-1,
        type="double",backingfile=""
      )
      bkname<-file.name(shared_mat)
    }else{
      fpath<-gsub("[^/]+$","",file_backed)
      fname<-gsub(fpath,"",file_backed)
      shared_mat<-big.matrix(
        nrow=nrow(mat)*2-1,ncol=nrow(mat)*2-1,type="double",
        backingfile=fname,backingpath=fpath,descriptorfile=paste0(fname,".desc")
      )
      bkname<-file_backed
    }

    cat("backing file: ",bkname,sep="",fill=T)
    edge<-greedyJoinEdgeBig(mat,weight=weight,pShared=shared_mat@address,verbose=verbose)
    rm(shared_mat)
    gc()
    if(file.exists(bkname)){
      unlink(bkname)
    }
    if(file.exists(paste0(bkname,".desc"))){
      unlink(paste0(bkname,".desc"))
    }
  }else{
    edge<-greedyJoinEdge(mat,weight=weight,verbose=verbose)
  }
  
  state_mrca<-charAncestral(edge[edge[,1]==2*nrow(mat)-2,1:ncol(mat)+3])
  if(collapse){
    edge<-collapseBranch(edge,ntip=nrow(mat))
    edge<-edge[edge[,3]!=-1,]
  }
  edge[,1:2]<-edge[,1:2]+1
  edge<-reOrderEdge(edge,nrow(mat)*2-1)
  phy$edge.length<-edge[,3]
  phy$edge<-edge[,1:2]
  temp<-unique(as.vector(phy$edge))
  temp<-temp[temp>nrow(mat)]
  phy$edge<-matrix(as.integer(factor(phy$edge,levels=c(1:nrow(mat),temp))),,2)
  phy$Nnode=max(phy$edge)-nrow(mat)
  phy$tip.label=rownames(mat)
  class(phy)<-"phylo"
  attr(phy,"order")<-"cladewise"
  state_ans<-rbind(state_mrca,edge[,1:ncol(mat)+3])
  state_ans<-state_ans[order(c(nrow(mat)+1,phy$edge[,2])),]
  
  return(list(phy=phy,ancestral=state_ans))
}

greedy.partial<-function(inmat,outmat=NULL,keep_cutoff=1,weight=NULL,noise=F,collapse=F,file_backed=NULL,verbose=F){
  phy<-greedy.wrap(inmat,weight=weight,noise=noise,collapse=collapse,file_backed=file_backed,verbose=verbose)
  if(is.null(outmat)){
    return(phy$phy)
  }
  node_mat<-matrix(phy$ancestral[(1:nrow(inmat))*(-1),],,ncol(inmat))
  phy<-phy$phy
  outmat_retain<-apply(outmat,1,function(x){
    mutationNested(x,node_mat[1,])
  })
  outmat_retain<-outmat_retain<keep_cutoff
  if(sum(outmat_retain)==0){
    return(list(phys=list(phy),category=rep(1,nrow(inmat)),isbroken=F))
  }else if(sum(outmat_retain)==1){
    outmat<-matrix(outmat[outmat_retain,],,ncol(inmat))
  }else{
    outmat<-outmat[outmat_retain,]
  }
  
  nodes_retain<-nodeKeep(edge=phy$edge,
                         Ntip=Ntip(phy),
                         node_mat=node_mat,
                         outmat=outmat,
                         keep_cutoff=keep_cutoff)
  nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
  broken_phys<-sub.multi.tree(phy,nodes_retain)
  
  category<-rep(0,nrow(inmat))
  i<-1
  while(i<=length(nodes_retain)){
    temp<-edge.subtree(phy,nodes_retain[i])
    category[temp[temp[,2]<=Ntip(phy),2]]<-i
    i<-i+1
  }
  
  names(category)<-rownames(inmat)
  return(list(phys=broken_phys,category=category,isbroken=T,origin_phy=phy,nodes_retain=nodes_retain))
}

#' phylogeny reconstruction using shared-mutation-join algorithm
#'
#' @param mat The mutation matrix.
#' @param bit_coded Whether the matrix is bit coded.
#' @param weight Weights of indels.
#' 
#' @examples
#' \dontrun{
#' tre<-shared.mutation.join(mat)
#'}
#'
#' @return a reconstructed phylogeny.

shared.mutation.join<-function(mat,bit_coded=F,weight=NULL){
  if(bit_coded){
    if(is.null(weight)){
      weight<-rep(1,ncol(mat)*32)
    }
  }else{
    mat<-regularize.mat(mat,mat)
    if(is.null(weight)){
      tmp<-mat2bit(mat,return.bit=F)
      weight<-matrix(1,max(diff(tmp$loci)),ncol(mat))
      rm(tmp)
    }
  }
  phy<-phy.join(inmat=t(mat),bit_coded=bit_coded,weight=weight)
  return(phy$phy)
}

h.join<-function(mat,splt,weight=NULL,collapse=F,verbose=F){
  mat<-regularize.mat(mat,mat)
  splt<-cbind(splt,0)
  k<-1
  output<-list()
  for(i in ncol(splt):2-1){
    for(j in setdiff(sort(unique(splt[,i])),0)){
      idx<-which(splt[,i]==j)
      submat<-gen.mat(mat[idx,],splt[idx,i+1])
      if(verbose){
        cat("Solving ",j,"th monogroup with ",nrow(submat$submat)," tips",sep="",fill=T)
      }
      output[[k]]<-greedy.partial(
        inmat=submat$submat,outmat=NULL,
        collapse=collapse,weight=weight,verbose=verbose
      )
      names(output)[k]<-paste0("M",j)
      k<-k+1
    }
  }
  return(output)
}

#' link clusters to a single phylogeny
#'
#' @param mg List of phylogenies.
#' 
#' @examples
#' \dontrun{
#' tre<-link.hcluster(phys)
#'}
#'
#' @return a phylogeny.

# link.hcluster<-function(mg){
#   output<-mg[["__M1"]]
#   while(length(grep("^__M",output$tip.label))!=0){
#     tip_idx<-grep("^__M",output$tip.label)[1]
#     output<-bind.tree(output,mg[[output$tip.label[tip_idx]]],where=tip_idx)
#   }
#   return(output)
# }

link.hcluster<-function(phys){
  tre<-phys[["__M1"]]
  while(length(grep("__M",tre$tip.label))>0){
    linkIdx<-grep("__M",tre$tip.label)
    edgeIdx<-match(linkIdx,tre$edge[,2])
    edgeIdx<-sort(edgeIdx)
    linkIdx<-tre$edge[edgeIdx,2]
    
    physAdd<-phys[tre$tip.label[linkIdx]]
    ntipAdd<-sapply(physAdd,Ntip)
    tipOffset<-diffinv(ntipAdd)+Ntip(tre)-length(linkIdx)
    nnodeAdd<-sapply(physAdd,Nnode)
    nodeOffset<-diffinv(nnodeAdd)[1:length(nnodeAdd)]
    lastNodeIdx<-sapply(edgeIdx,function(i){
      max(tre$edge[1:i,1])
    })
    nodeOffset<-nodeOffset+lastNodeIdx+sum(ntipAdd)-length(ntipAdd)
    tipTrans<-sapply(1:Ntip(tre),function(i){
      i-sum(i>linkIdx)
    })
    nodeTrans<-sapply(1:Nnode(tre)+Ntip(tre),function(i){
      i+sum(ntipAdd)-length(ntipAdd)+sum(nnodeAdd[i>lastNodeIdx])
    })
    trans<-c(tipTrans,nodeTrans)
    
    edgeIdx<-c(0,edgeIdx)
    output<-list()
    output$edge<-matrix(0,2,Ntip(tre)+Nnode(tre)+sum(ntipAdd)-length(ntipAdd)+sum(nnodeAdd)-1)
    output$edge.length<-rep(0,Ntip(tre)+Nnode(tre)+sum(ntipAdd)-length(ntipAdd)+sum(nnodeAdd)-1)
    offset<-0
    
    for(i in 1:length(physAdd)){
      ntip<-Ntip(physAdd[[i]])
      nnode<-Nnode(physAdd[[i]])
      edgeTmp<-physAdd[[i]]$edge
      edgeTmp[edgeTmp>ntip]<-edgeTmp[edgeTmp>ntip]-ntip+nodeOffset[i]
      edgeTmp[edgeTmp<=ntip]<-edgeTmp[edgeTmp<=ntip]+tipOffset[i]
      
      len<-edgeIdx[i+1]-edgeIdx[i]
      output$edge[,1:len+offset]<-trans[t(tre$edge[1:len+edgeIdx[i],])]
      output$edge[2,len+offset]<-edgeTmp[1,1]
      if(!is.null(tre$edge.length)){
        output$edge.length[1:len+offset]<-tre$edge.length[1:len+edgeIdx[i]]
      }
      offset<-offset+len
      output$edge[,1:nrow(edgeTmp)+offset]<-t(edgeTmp)
      if(!is.null(physAdd[[i]]$edge.length)){
        output$edge.length[1:nrow(edgeTmp)+offset]<-physAdd[[i]]$edge.length
      }
      offset<-offset+nrow(edgeTmp)
    }
    if(offset<ncol(output$edge)){
      output$edge[,(offset+1):ncol(output$edge)]<-trans[t(tre$edge[(edgeIdx[i+1]+1):nrow(tre$edge),])]
      if(!is.null(tre$edge.length)){
        output$edge.length[(offset+1):ncol(output$edge)]<-tre$edge.length[(edgeIdx[i+1]+1):nrow(tre$edge)]
      }
    }
    output$tip.label<-c(tre$tip.label[setdiff(1:Ntip(tre),linkIdx)],
                        unname(unlist(lapply(physAdd,function(x) x$tip.label))))
    output$Nnode<-tre$Nnode+sum(unlist(lapply(physAdd,Nnode)))
    output$edge<-t(output$edge)
    class(output)<-"phylo"
    attr(output,"order")<-"cladewise"
    tre<-output    
  }
  return(tre)
}

edge.subtree<-function(phy,node){
  if(node==length(phy$tip.label)+1){
    return(cbind(phy$edge,phy$edge.length))
  }
  edge<-cbind(phy$edge,phy$edge.length)
  edge<-edge[which(edge[,1]==node)[1]:nrow(edge),]
  if(min(edge[,1])<node){
    edge<-edge[1:(which(edge[,1]<node)[1]-1),]
  }
  return(edge)
}

edge.wrap<-function(edge,phy){
  output<-list()
  if(ncol(edge)>2){
    output$edge.length<-edge[,3]
  }
  edge<-edge[,1:2]
  output$tip.label<-phy$tip.label[sort(edge[edge[,2]<=length(phy$tip.label),2])]
  edge<-matrix(as.integer(factor(edge)),,2)
  output$Nnode<-max(edge)-length(output$tip.label)
  output$edge<-edge
  class(output)<-"phylo"
  attr(output,"order")<-"cladewise"
  return(output)
}

sub.one.tree<-function(phy,node){
  edge<-edge.subtree(phy,node)
  return(edge.wrap(edge,phy))
}

sub.multi.tree<-function(phy,nodes){
  output<-lapply(nodes,function(i){
    edge.wrap(edge.subtree(phy,i),phy)
  })
  return(output)
}

nodes.break<-function(phy,nodes_state){
  output<-list()
  k<-1
  exist_subtree<-sum(nodes_state)
  while(exist_subtree){
    node_idx<-which(nodes_state)[1]+length(phy$tip.label)
    edge<-edge.subtree(phy,node_idx)
    nodes_state[unique(edge[,1])-length(phy$tip.label)]<-F
    output[[k]]<-edge.wrap(edge,phy)
    k<-k+1
    exist_subtree<-sum(nodes_state)
  }
  return(output)
}

gen.mat<-function(mat,splt,bit_coded=F,return.mat=T){
  mg_idx<-setdiff(unique(splt),0)
  if(length(mg_idx)==0){
    submat_idx<-1:nrow(mat)
    if(return.mat){
      return(list(submat=mat,idx=submat_idx,splt=rep(0,nrow(mat))))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=rep(0,nrow(mat))))
    }
  }else{
    if(return.mat){
      if(bit_coded){
        output<-rbind(
          output<-mat[splt==0,],
          t(sapply(mg_idx,function(x){
            bitAncestral(matrix(mat[splt==x,],,ncol(mat)))
          }))
        )
      }else{
        output<-rbind(
          output<-mat[splt==0,],
          t(sapply(mg_idx,function(x){
            charAncestral(matrix(mat[splt==x,],,ncol(mat)))
          }))
        )
      }
      rownames(output)<-c(rownames(mat)[splt==0],paste0("__M",mg_idx))
    }
    offset<-sum(splt==0)
    submat_idx<-sapply(splt,function(x){
      if(x==0) return(0)
      else return(match(x,mg_idx)+offset)
    })
    submat_idx[submat_idx==0]<-1:offset
    if(return.mat){
      return(list(submat=output,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }
  }
}

gen.mat2<-function(mat,ansMat=NULL,splt,bit_coded=F,return.mat=T){
  ## col:tip row:feature
  mg_idx<-setdiff(unique(splt),0)
  if(length(mg_idx)==0){
    submat_idx<-1:ncol(mat)
    if(return.mat){
      return(list(submat=mat,idx=submat_idx,splt=rep(0,ncol(mat))))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=rep(0,ncol(mat))))
    }
  }else{
    if(return.mat){
      output<-cbind(
        output<-mat[,splt==0],
        ansMat[,paste0("__M",mg_idx)]
      )
      colnames(output)<-c(colnames(mat)[splt==0],paste0("__M",mg_idx))
    }
    offset<-sum(splt==0)
    submat_idx<-sapply(splt,function(x){
      if(x==0) return(0)
      else return(match(x,mg_idx)+offset)
    })
    submat_idx[submat_idx==0]<-1:offset
    if(return.mat){
      return(list(submat=output,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }
  }
}

gen.mat3<-function(mat,inIdx,ansMat=NULL,splt,return.mat=T){
  ##row:tip col:feature
  mg_idx<-setdiff(unique(splt),0)
  if(length(mg_idx)==0){
    submat_idx<-1:length(inIdx)
    if(return.mat){
      return(list(submat=mat[inIdx,],idx=submat_idx,splt=rep(0,length(inIdx))))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=rep(0,length(inIdx))))
    }
  }else{
    if(return.mat){
      output<-rbind(
        output<-mat[inIdx[splt==0],],
        ansMat[paste0("__M",mg_idx),]
      )
      rownames(output)<-c(rownames(mat)[inIdx][splt==0],paste0("__M",mg_idx))
    }
    offset<-sum(splt==0)
    submat_idx<-sapply(splt,function(x){
      if(x==0) return(0)
      else return(match(x,mg_idx)+offset)
    })
    submat_idx[submat_idx==0]<-1:offset
    if(return.mat){
      return(list(submat=output,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }
  }
}


#' reconstruct the phylogenies using partial join algorithm
#'
#' @param mat The mutation matrix.
#' @param splt Split of data from `split.data`.
#' @param keep_cutoff The cutoff to break problematic clades.
#' @param weight Weights of indels.
#' @param bit_coded Whether the matrix is bit coded.
#' @param prefix The prefix of temporary file.
#' @param cuda_path The path of CUDA script. If set to `NULL`, do not use CUDA.
#' @param cuda_threshold Cluster larger than this value will be accelerated by CUDA.
#' @param ncore The number of cores used.
#' @param verbose Whether print detailed information.
#' 
#' @examples
#' \dontrun{
#' phys<-h.partial.join(imputed,splt,noise=F,verbose=F)
#'}
#'
#' @return a list of phylogenies.

h.partial.join<-function(
    mat,splt,bit_coded=F,keep_cutoff=1,weight=NULL,
    prefix="./",cuda_path=NULL,cuda_threshold=500,
    ncore=1,verbose=T
){
  ancestralType=1
  sizeOfTips=rep(1,nrow(mat))
  minSize=0
  minRatio=2
  
  partial.join<-function(j){
    inIdx<-which(splt[,layer]==j)
    outIdx<-which(splt[,layer]!=j)
    inmat<-gen.mat3(mat,inIdx,ansMat,splt[inIdx,layer+1])
    tip2pseudo<-inmat$idx
    if(!bit_coded){
      tipSize<-sapply(1:nrow(inmat$submat),function(x){
        sum(sizeOfTips[inIdx][tip2pseudo==x])
      })
      tipSize<-as.integer(tipSize)
      # tipSize<-calPseudoSize(sizeOfTips[inIdx],tip2pseudo,nrow(inmat$submat))
    }
    inmat<-t(inmat$submat)
    tipLabel<-colnames(inmat)
    
    ##using CUDA to join or not
    if(ncol(inmat)<cuda_threshold||is.null(cuda_path)){
      phyData<-phy.join(
        inmat=inmat,
        bit_coded=bit_coded,
        weight=w,
        ancestralType=ancestralType,
        tipSize=tipSize,
        minSize=minSize,
        minRatio=minRatio
      )
      rm(inmat)
      phy<-phyData$phy
      node_mat<-phyData$node_mat
      rm(phyData)
      if(length(outIdx)==0){
        return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
      }
      
      tmpSplt<-splt[outIdx,layer+1]
      idx1<-outIdx[tmpSplt==0]
      mgIdx<-setdiff(unique(tmpSplt),0)
      if(length(mgIdx)>0){
        idx2<-match(paste0("__M",setdiff(unique(tmpSplt),0)),rownames(ansMat))
      }else{
        idx2<-integer()
      }
      if(bit_coded){
        idx1f<-outmatFilterBit2(mat,idx1-1L,node_mat[,1],cutoff=keep_cutoff)+1L
        idx2f<-outmatFilterBit2(ansMat,idx2-1L,node_mat[,1],cutoff=keep_cutoff)+1L
      }else{
        idx1f<-outmatFilterIndel2(mat,idx1-1L,node_mat[,1],cutoff=keep_cutoff)+1L
        idx2f<-outmatFilterIndel2(ansMat,idx2-1L,node_mat[,1],cutoff=keep_cutoff)+1L
      }
      if(length(idx1f)+length(idx2f)==0){
        return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
      }else{
        outmat<-t(rbind(mat[idx1f,],ansMat[idx2f,]))
      }
      
      if(bit_coded){
        ## the results are copied from edge, the idx of thy$edge is 1-based, no need to plus one
        nodes_retain<-nodeKeepBit2(
          edge=phy$edge,
          Ntip=Ntip(phy),
          node_mat=node_mat,
          outmat=outmat,
          keep_cutoff=keep_cutoff
        )
      }else{
        nodes_retain<-nodeKeepIndel(
          edge=phy$edge,
          Ntip=Ntip(phy),
          node_mat=node_mat,
          outmat=outmat,
          keep_cutoff=keep_cutoff
        )
      }
      rm(outmat)
      nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
      node_mat<-t(node_mat[,nodes_retain-Ntip(phy)])
      gc()
    }else{
      if(bit_coded){
        write.mat2bin(inmat,file=paste0(prefix,"__M",j,".tbmat",sep=""))
        rm(inmat)
        system(paste0(cuda_path,"joinEdgeBitReverse ",prefix,"__M",j," ",wfile," ",ancestralType))
      }else{
        write.mat2bin(inmat+1L,file=paste0(prefix,"__M",j,".tbmat",sep=""))
        write.mat2bin(matrix(tipSize,,1),file=paste0(prefix,"__M",j,".tipSize",sep=""))
        rm(inmat)
        system(paste0(cuda_path,"joinEdgeIndel ",prefix,"__M",j," ",wfile," ",nrow(w)," ",minSize," ",minRatio))
      }
      unlink(paste0(prefix,"__M",j,".tbmat"))
      unlink(paste0(prefix,"__M",j,".tipSize",sep=""))
      edge<-read.bin2mat(paste0(prefix,"__M",j,".edge2"))
      phy<-edge2phy(edge,tipLabel)
      node_mat<-read.bin2mat(paste0(prefix,"__M",j,".nodeState2",sep=""))
      if(!bit_coded){
        node_mat<-node_mat-1L
      }
      if(length(outIdx)==0){
        unlink(paste0(prefix,"__M",j,".edge2"))
        unlink(paste0(prefix,"__M",j,".nodeState2"))
        return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
      }
      
      tmpSplt<-splt[outIdx,layer+1]
      idx1<-outIdx[tmpSplt==0]
      mgIdx<-setdiff(unique(tmpSplt),0)
      if(length(mgIdx)>0){
        idx2<-match(paste0("__M",setdiff(unique(tmpSplt),0)),rownames(ansMat))
      }else{
        idx2<-integer()
      }
      if(bit_coded){
        idx1f<-outmatFilterBit2(mat,idx1-1L,node_mat[,1],cutoff=keep_cutoff)+1L
        idx2f<-outmatFilterBit2(ansMat,idx2-1L,node_mat[,1],cutoff=keep_cutoff)+1L
      }else{
        idx1f<-outmatFilterIndel2(mat,idx1-1L,node_mat[,1],cutoff=keep_cutoff)+1L
        idx2f<-outmatFilterIndel2(ansMat,idx2-1L,node_mat[,1],cutoff=keep_cutoff)+1L
      }
      if(length(idx1f)+length(idx2f)==0){
        unlink(paste0(prefix,"__M",j,".edge2"))
        unlink(paste0(prefix,"__M",j,".nodeState2"))
        return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
      }else{
        outmat<-t(rbind(mat[idx1f,],ansMat[idx2f,]))
      }
      
      ########using CUDA to break phy or not
      if(ncol(outmat)<128){
        unlink(paste0(prefix,"__M",j,".edge2"))
        unlink(paste0(prefix,"__M",j,".nodeState2"))
        if(bit_coded){
          ## the results are copied from edge, the idx of thy$edge is 1-based, no need to plus one
          nodes_retain<-nodeKeepBit2(
            edge=phy$edge,
            Ntip=Ntip(phy),
            node_mat=node_mat,
            outmat=outmat,
            keep_cutoff=keep_cutoff
          )
        }else{
          nodes_retain<-nodeKeepIndel(
            edge=phy$edge,
            Ntip=Ntip(phy),
            node_mat=node_mat,
            outmat=outmat,
            keep_cutoff=keep_cutoff
          )
        }
        rm(outmat)
        nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
        node_mat<-t(node_mat[,nodes_retain-Ntip(phy)])
        gc()
      }else{
        if(bit_coded){
          write.mat2bin(outmat,file=paste0(prefix,"__M",j,".tboutmat",sep=""))
          rm(outmat)
          system(paste0(cuda_path,"nodeKeepBit ",prefix,"__M",j," ",keep_cutoff))
        }else{
          write.mat2bin(outmat+1L,file=paste0(prefix,"__M",j,".tboutmat",sep=""))
          rm(outmat)
          system(paste0(cuda_path,"nodeKeepIndel ",prefix,"__M",j," ",keep_cutoff))
        }
        if(file.size(paste0(prefix,"__M",j,".retainNode"))==0){
          nodes_retain<-integer()
        }else{
          nodes_retain<-read.table(file=paste0(prefix,"__M",j,".retainNode"))[,1]+1
          nodes_retain<-sort(nodes_retain)
        }
        unlink(paste0(prefix,"__M",j,".edge2"))
        unlink(paste0(prefix,"__M",j,".nodeState2"))
        unlink(paste0(prefix,"__M",j,".tboutmat"))
        unlink(paste0(prefix,"__M",j,".retainNode"))
        node_mat<-t(node_mat[,nodes_retain-Ntip(phy)])
        gc()
      }
    }
    if(length(nodes_retain)==0){
      return(list(phys=list(),category=rep(0,length(inIdx)),ansState=node_mat,isbroken=T,inIdx=inIdx))
    }else{
      tmp<-phy.break(phy,nodes_retain)
      broken_phys<-tmp$phys
      category<-tmp$category
      return(list(phys=broken_phys,category=category[tip2pseudo],ansState=node_mat,isbroken=T,inIdx=inIdx))
    }
  }
  wfile<-paste0(prefix,".w",sep="")
  if(!bit_coded){
    mat<-regularize.mat(mat,mat)
    if(is.null(weight)){
      tmp<-mat2bit(mat,return.bit=F)
      w<-matrix(1,max(diff(tmp$loci)),ncol(mat))
      rm(tmp)
    }else{
      w<-weight
    }
    writeBin(as.numeric(t(w)),con=wfile)
  }else{
    if(is.null(weight)){
      w<-rep(1,ncol(mat)*32)
    }else{
      w<-weight
    }
    writeBin(as.numeric(w),con=wfile)
  }
  
  splt<-cbind(splt,0)
  max_mg_idx<-max(splt)
  phys<-list()
  ansMat<-matrix(integer(),0,ncol(mat))
  for(layer in ncol(splt):2-1){
    clustIdx<-setdiff(sort(unique(splt[,layer])),0)
    set.seed(1)
    if(length(clustIdx)>1){
      clustIdx<-sample(clustIdx,length(clustIdx))
    }
    cat(as.character(Sys.time()),"solving layer ",layer,fill=T)
    ncoreRun<-min(ncore,length(clustIdx))
    if(ncoreRun>1){
      cl<-makeForkCluster(ncoreRun)
      partial_res_collected<-parLapply(cl,clustIdx,partial.join)
      stopCluster(cl)
    }else{
      partial_res_collected<-list()
      for(j in clustIdx){
        partial_res_collected<-c(partial_res_collected,list(partial.join(j)))
      }
    }
    cat(as.character(Sys.time()),"summerizing layer ",layer,fill=T)
    nPhys<-sapply(partial_res_collected,function(x){
      if(x$isbroken){
        return(length(x$phys))
      }else{
        return(0)
      }
    })
    offset<-diffinv(nPhys)+max_mg_idx
    tmp<-unlist(lapply(1:length(clustIdx),function(i){
      if(partial_res_collected[[i]]$isbroken){
        if(nPhys[i]==0){
          return(integer())
        }else{
          return(1:nPhys[i]+offset[i])
        }
      }else{
        return(clustIdx[i])
      }
    }))
    newPhys<-unlist(lapply(partial_res_collected,function(x) x$phys),recursive=F)
    if(length(tmp)>0){
      names(newPhys)<-paste0("__M",tmp)
    }
    phys<-c(phys,newPhys)
    
    newAnsMat<-lapply(partial_res_collected,function(x) t(x$ansState))
    newAnsMat<-t(matrix(unlist(newAnsMat),ncol(mat),))
    if(length(tmp)>0){
      rownames(newAnsMat)<-paste0("__M",tmp)
    }
    ansMat<-rbind(ansMat,newAnsMat)
    
    replaceIdx<-unlist(lapply(partial_res_collected,function(x){
      if(x$isbroken){
        return(x$inIdx)
      }else{
        return(integer())
      }
    }))
    new_group<-unlist(lapply(1:length(clustIdx),function(i){
      if(partial_res_collected[[i]]$isbroken){
        tmp<-partial_res_collected[[i]]$category
        tmp[tmp!=0]<-tmp[tmp!=0]+offset[i]
        return(tmp)
      }else{
        return(integer())
      }
    }))
    temp<-splt[replaceIdx,layer+1]
    temp[new_group>0]<-new_group[new_group>0]
    splt[replaceIdx,layer]<-temp
    max_mg_idx<-max(offset)
    cat(as.character(Sys.time())," layer ",layer," finished",fill=T)
  }
  unlink(wfile)
  return(phys)
}

#' calculate the number of evolution events of the phylogeny.
#'
#' @param phy The phylogeny object.
#' @param mat The mutation matrix.
#' @param bit_coded Whether the matrix is bit coded.
#' 
#' @examples
#' \dontrun{
#' nEvo<-cal.nevo(tre,mat)
#'}
#'
#' @return an bit-coded matrix marking the evolution events in the nodes of the phylogeny.

cal.nevo<-function(phy,mat,bit_coded=F){
  if(bit_coded){
    nEvo<-ancestralBitReverse(t(mat),phy$edge,Ntip(phy),Nnode(phy))
  }else{
    nEvo<-charAncestralEC(t(mat),phy$edge,tipSize=rep(1,nrow(mat)),Nnode(phy),minSize=0,minRatio=0)
  }
  return(nEvo)
}

add.nEvo<-function(phy,nEvo,weight=NULL){
  if(is.null(weight)){
    weight<-rep(1,nrow(nEvo)*32)
  }
  phy$edge.length<-apply(nEvo,2,locusSUM,weight=weight)[phy$edge[,2]]
  return(phy)
}

#' add the length of edges to the phylogeny and collapse nodes without evolution.
#'
#' @param phy The phylogeny object.
#' @param mat The mutation matrix.
#' @param bit_coded Whether the matrix is bit coded.
#' @param collapse Whether to collapse the nodes. If set to False, the function just add the length of edges to the phylogeny.
#' 
#' @examples
#' \dontrun{
#' tre<-collapse.infoless(tre,mat)
#'}
#'
#' @return a phylogeny.

collapse.infoless<-function(phy,mat,bit_coded=F,collapse=T){
  if(bit_coded){
    nEvo<-ancestralBitReverse(t(mat),phy$edge,Ntip(phy),Nnode(phy))
  }else{
    nEvo<-charAncestralEC(t(mat),phy$edge,rep(1,Ntip(phy)),Nnode(phy),minSize=0,minRatio=0)
  }
  phy$edge.length<-nEvo$edgeLen
  if(collapse){
    phy<-di2multi(phy)
  }
  return(phy)
}

plot.nevo<-function(mat,nEvo1,nEvo2=NULL,cols,weight=NULL){
  binary<-mat2bit(mat,return.bit=F)
  nloci<-rev(binary$loci)[1]
  bin_size<-diff(binary$loci)
  if(!is.null(weight)){
    weight<-unlist(lapply(which(bin_size!=0),function(i){
      weight[1:bin_size[i],i]
    }))
  }else{
    weight<-rep(1,nloci)
  }
  xx<-integer()
  yy<-integer()
  cols_idx<-integer()
  for(i in 2:length(binary$loci)-1){
    len<-binary$loci[i+1]-binary$loci[i]
    if(len!=0){
      yy<-c(yy,(1:len-1)*(1))
      xx<-c(xx,rep(i-1,len))
      cols_idx<-c(cols_idx,1:len+2)
    }
  }
  for(i in 1:length(bin_size)){
    if(bin_size[i]!=0){
      yy<-c(yy,(1:bin_size[i]-1)*(1))
      xx<-c(xx,rep(i-1,bin_size[i]))
      cols_idx<-c(cols_idx,1:bin_size[i]+2)
    }
  }
  if(is.null(nEvo2)){
    total_nEvo<-sum(apply(nEvo1,2,locusSUM,weight=weight))
  }else{
    total_nEvo<-c(
      sum(apply(nEvo2,2,locusSUM,weight=weight)),
      sum(apply(nEvo1,2,locusSUM,weight=weight))
    )
  }
  plot(x=xx,y=yy,col=cols[cols_idx],pch=15,
       xaxt="n",yaxt="n",main=paste0("No. evolution: ",round(total_nEvo,3)),
       xlab="",ylab="indels")
  axis(side=1,labels=paste0("locus ",1:ncol(mat)),at=1:ncol(mat)-1,las=2,cex.axis=0.5)
  axis(side=2,labels=1:max(diff(binary$loci)),at=1:max(diff(binary$loci))-1,las=1,cex.axis=0.5)
  if(is.null(nEvo2)){
    text(x=xx,y=yy,labels=rowSums(apply(nEvo1,2,int2binary,output_len=nloci)),cex=0.5)
  }else{
    text(x=xx,y=yy,
         labels=rowSums(apply(nEvo2,2,int2binary,output_len=nloci))-rowSums(apply(nEvo1,2,int2binary,output_len=nloci)),
         cex=0.5)
  }
}

cal.mono<-function(phy,return.full=F){
  ntip<-length(phy$tip.label)
  output<-matrix(F,ntip,ntip)
  diag(output)<-T
  output<-rbind(output,matrix(F,phy$Nnode,ntip))
  for(i in nrow(phy$edge):1){
    output[phy$edge[i,1],]<-output[phy$edge[i,1],]|output[phy$edge[i,2],]
  }
  if(return.full){
    return(output)
  }else{
    return(output[1:phy$Nnode+ntip,])
  }
}

summerize.evo<-function(phy,nEvo){
  output<-matrix(0,length(phy$tip.label),nrow(nEvo))
  mono_mat<-cal.mono(phy,return.full=T)
  apply(nEvo,1,function(x){
    if(sum(x)==1)
      return(mono_mat[which(x),])
    else
      colSums(mono_mat[which(x),]*(1:sum(x)))
  })
}

#' insert duplicated OTUs back to the phylogeny
#'
#' @param phy The phylogeny object.
#' @param dup_idx The indices of duplicated OTUs from `find.dup`.
#' @param label Tip label of the full phylogeny.
#' 
#' @examples
#' \dontrun{
#' tre<-insert.dup(tre,dup_idx)
#'}
#'
#' @return a phylogeny.

insert.dup<-function(phy,dup_idx,label=paste0("C",1:length(dup_idx))){
  phy$node.label<-NULL
  retain<-which(dup_idx==1:length(dup_idx))
  removed<-which(dup_idx!=1:length(dup_idx))
  ndup<-length(removed)
  anchor<-sort(unique(dup_idx[removed]))
  edge<-phy$edge
  edge[edge>Ntip(phy)]<-edge[edge>Ntip(phy)]+ndup
  edge[edge<=Ntip(phy)]<-retain[edge[edge<=Ntip(phy)]]
  node_offset<-Ntip(phy)+Nnode(phy)+ndup
  e1<-unlist(apply(edge,1,function(x){
    if(x[2]%in%anchor){
      return(c(x[1],rep(node_offset+match(x[2],anchor),sum(dup_idx==x[2]))))
    }else{
      return(x[1])
    }
  }))
  e2<-unlist(apply(edge,1,function(x){
    if(x[2]%in%anchor){
      return(c(node_offset+match(x[2],anchor),which(x[2]==dup_idx)))
    }else{
      return(x[2])
    }
  }))
  if(!is.null(phy$edge.length)){
    phy$edge.length<-unlist(sapply(1:nrow(edge),function(i){
      if(edge[i,2]%in%anchor){
        return(c(phy$edge.length[i],rep(0,sum(dup_idx==edge[i,2]))))
      }else{
        return(phy$edge.length[i])
      }
    }))
  }
  edge<-unname(cbind(e1,e2))
  edge[edge>Ntip(phy)+ndup]<-match(edge[edge>Ntip(phy)+ndup],unique(edge[edge>Ntip(phy)+ndup]))+Ntip(phy)+ndup
  
  phy$edge<-edge
  phy$tip.label<-label
  phy$Nnode<-phy$Nnode+length(anchor)
  return(phy)
}

#' re-order the tip of the phylogeny.
#'
#' @param phy The phylogeny object.
#' @param odr a character vector specifies the order of tips
#' 
#' @examples
#' \dontrun{
#' tre<-order.tip(tre,rev(tre$tip.label))
#'}
#'
#' @return a phylogeny.

order.tip<-function(phy,odr){
  phy<-collapse.singles(phy)
  ntip<-Ntip(phy)
  phy$edge[phy$edge[,2]<=ntip,2]<-match(phy$tip.label[phy$edge[phy$edge[,2]<=ntip,2]],odr)
  phy$tip.label<-odr
  return(phy)
}

#' calculate the number of OTUs after splitting the data.
#'
#' @param splt The results from `split.data`.
#' 
#' @examples
#' \dontrun{
#' head(splt.summerize(splt),n=20)
#'}
#'
#' @return an dataframe containing the layer, the cluster index and the number of OTUs.

splt.summerize<-function(splt){
  splt<-cbind(splt,0)
  output<-matrix(,0,3)
  for(i in 2:ncol(splt)-1){
    for(j in setdiff(sort(unique(splt[,i])),0)){
      ntip<-length(setdiff(unique(splt[splt[,i]==j,i+1]),0))
      ntip<-ntip+sum(splt[splt[,i]==j,i+1]==0)
      output<-rbind(output,c(i,j,ntip))
    }
  }
  output<-as.data.frame(output)
  colnames(output)<-c("Layer","Monogrup","No. tips")
  rownames(output)<-NULL
  return(output)
}

gen.weight<-function(mat,nEvo){
  binary<-mat2bit(mat,return.bit=F)
  temp<-diff(binary$loci)
  nEvo<-rowSums(apply(nEvo,2,int2binary))
  weight<-sapply(1:ncol(mat),function(i){
    output<-rep(0,max(temp))
    if(temp[i]!=0){
      output[1:temp[i]]<-1/nEvo[1:temp[i]+binary$loci[i]]
    }
    return(output)
  })
  weight[weight!=0]<-weight[weight!=0]/mean(weight[weight!=0])
  return(weight)
}

random.weight<-function(mat,scale=0.01){
  binary<-mat2bit(mat,return.bit=F)
  temp<-diff(binary$loci)
  weight<-matrix((runif(max(temp)*ncol(mat))*2-1)*scale+1,,ncol(mat))
  return(weight)
}

node.ident<-function(phy,ident){
  ntip<-Ntip(phy)
  node_member<-nodeDes(phy$edge,Nnode(phy),Ntip(phy))
  node_ident<-apply(node_member,2,function(x){
    tmp<-ident[int2binary(x,output_len=ntip)]
    if(length(unique(tmp))==1){
      return(unique(tmp))
    }else{
      return(0)
    }
  })
  return(node_ident)
}

find.des<-function(phys,id){
  output<-phys[[id]]$tip.label
  while(length(grep("M",output))!=0){
    temp<-grep("M",output,value=T)
    temp<-sapply(phys[temp],function(x){
      x$tip.label
    })
    output<-c(grep("C",output,value=T),unlist(temp))
  }
  return(unname(output))
}

join.sample.tree<-function(multi_tree,tag,ident,sample_idx){
  batch<-unique(ident[sample_idx])
  phy<-list()
  phy$tip.label<-paste0("M",batch)
  phy$edge<-cbind(length(batch)+1,1:length(batch))
  phy$Nnode<-1
  attr(phy,"class")="phylo"
  attr(phy,"order")="cladewise"
  for(i in 1:length(batch)){
    temp<-multi_tree[[batch[i]]]
    temp$tip.label<-tag[ident==batch[i]]
    temp<-keep.tip(temp,which(which(ident==batch[i])%in%sample_idx))
    tip_idx<-which(phy$tip.label==paste0("M",batch[i]))
    phy<-bind.tree(phy,temp,where=tip_idx)
  }
  return(phy)
}

#' find duplicated OTUs in the mutation matrix.
#'
#' @param mat The mutation matrix.
#' @param ncore The number of threads used.
#' 
#' @examples
#' \dontrun{
#' dup_idx<-find.dup(mat,ncore=10)
#'}
#'
#' @return a integer vector. The value of ith element is the index of OTU shared the same mutation profiles with the ith OTU.

find.dup<-function(mat,bit_coded=F,ncore=1){
  if(!bit_coded){
    binary<-mat2bit(mat,code_unmutated=T)$binary
    if(ncore>1){
      cl<-makeForkCluster(ncore)
      output<-parSapply(cl,1:nrow(binary),function(i){
        markDup(i-1,binary)
      })
      output<-parSapply(cl,1:nrow(binary),function(i){
        transitionRoot(output,i-1)
      })
      stopCluster(cl)
    }else{
      output<-sapply(1:nrow(binary),function(i){
        markDup(i-1,binary)
      })
      output<-sapply(1:nrow(binary),function(i){
        transitionRoot(output,i-1)
      })
    }
  }else{
    if(ncore==1){
      output<-sapply(1:nrow(mat),function(i){
        markDupBit(i-1,mat)
      })
    }else{
      idx<-sample(1:nrow(mat),nrow(mat))
      cl<-makeForkCluster(ncore)
      output<-parSapply(cl,idx,function(i){
        markDupBit(i-1,mat)
      })
      stopCluster(cl)
      output<-output[order(idx)]
    }
  }
  return(output+1)
}

combine.dup<-function(dup_idx1,dup_idx2){
  retain<-which(dup_idx1==1:length(dup_idx1))
  dup_idx1<-retain[dup_idx2[match(dup_idx1,retain)]]
  return(dup_idx1)
}

#' calculate the bootstrap values of the phylogeny.
#'
#' @param phy The phylogeny object.
#' @param bsphy List of bootstrap trees.
#' @param ncore The number of threads used.
#' 
#' @examples
#' \dontrun{
#' bs_value<-cal.bs(tre,bstrees,10)
#'}
#'
#' @return An integer vector containing the bootstrap value of nodes in the phylogeny.

cal.bs<-function(phy,bsphy,ncore=1){
  if(ncore==1){
    bs<-prop.clades(phy,bsphy)
    bs[is.na(bs)]<-0
    return(bs[-1])
  }else{
    cl<-makeForkCluster(ncore)
    bs<-parSapply(cl,bsphy,function(x){
      prop.clades(phy,x)
    })
    stopCluster(cl)
    return(rowSums(bs,na.rm=T)[-1])
  }

}

read.phy<-function(file){
  mat<-read.table(file=file,colClass="character")[c(-1,-2),]
  output<-t(sapply(mat[,2],function(x){
    unname(unlist(strsplit(x,split="")))
  }))
  attr(output,"dimnames")<-NULL
  output<-matrix(as.integer(output),nrow(output),ncol(output))
  rownames(output)<-mat[,1]
  return(output)
}

#' turn an integer matrix to a bit-coded matrix
#'
#' @param mat The mutation matrix containing only 0 and 1.
#' @param check Whether to remove columns without information.
#' 
#' @examples
#' \dontrun{
#' mat<-to.bit(mat)
#'}
#'
#' @return an bit-coded matrix.

to.bit<-function(mat,check=F){
  if(check){
    col_idx<-which(apply(mat,2,function(x) sum(x==1)>1))
    mat<-mat[,col_idx]
  }
  rname<-rownames(mat)
  nloci<-ncol(mat)
  mat<-matCoding(mat,0:nloci)
  attr(mat,"bit_coded")<-T
  attr(mat,"nloci")<-nloci
  rownames(mat)<-rname
  return(mat)
}

regularize.mat.bit<-function(mat){
  nloci<-attr(mat,"nloci")
  return(t(apply(mat,1,int2binary)))
}

node.ancestral.bit<-function(phy,mat){
  output<-matrix(integer(),Ntip(phy)+Nnode(phy),ncol(mat))
  output[1:length(phy$tip.label),]<-mat
  for(i in Nnode(phy):1+Ntip(phy)){
    idx<-phy$edge[phy$edge[,1]==i,2]
    output[i,]<-bitAncestral(output[idx,])
  }
  output<-output[1:Nnode(phy)+Ntip(phy),]
  attr(output,"bit_coded")<-T
  attr(output,"nloci")<-attr(mat,"nloci")
  return(output)
}

greedy.wrap.bit<-function(mat,weight=NULL,noise=F,file_backed=NULL,verbose=F){
  if(is.null(weight)){
    weight<-rep(1,ncol(mat)*32)
  }
  if(noise){
    weight<-weight+runif(length(weight))/1e4
  }
  phy<-list()
  if(!is.null(file_backed)){
    if(file_backed==""){
      shared_mat<-big.matrix(
        nrow=nrow(mat)*2-1,ncol=nrow(mat)*2-1,
        type="double",backingfile=""
      )
      bkname<-file.name(shared_mat)
    }else{
      fpath<-gsub("[^/]+$","",file_backed)
      fname<-gsub(fpath,"",file_backed)
      shared_mat<-big.matrix(
        nrow=nrow(mat)*2-1,ncol=nrow(mat)*2-1,type="double",
        backingfile=fname,backingpath=fpath,descriptorfile=paste0(fname,".desc")
      )
      bkname<-file_backed
    }
    cat("backing file: ",bkname,sep="",fill=T)
    edge<-greedyJoinEdgeBitBig(mat,weight=weight,pShared=shared_mat@address,verbose=verbose)
    rm(shared_mat)
    gc()
    if(file.exists(bkname)){
      unlink(bkname)
    }
    if(file.exists(paste0(bkname,".desc"))){
      unlink(paste0(bkname,".desc"))
    }
  }else{
    edge<-greedyJoinEdgeBit(mat,weight=weight,verbose=verbose)
  }
  
  state_mrca<-bitAncestral(edge[edge[,1]==2*nrow(mat)-2,1:ncol(mat)+3])
  edge[,1:2]<-edge[,1:2]+1
  edge<-reOrderEdge(edge,nrow(mat)*2-1)
  phy$edge.length<-edge[,3]
  phy$edge<-edge[,1:2]
  temp<-unique(as.vector(phy$edge))
  temp<-temp[temp>nrow(mat)]
  phy$edge<-matrix(as.integer(factor(phy$edge,levels=c(1:nrow(mat),temp))),,2)
  phy$Nnode=max(phy$edge)-nrow(mat)
  phy$tip.label=rownames(mat)
  class(phy)<-"phylo"
  attr(phy,"order")<-"cladewise"
  state_ans<-rbind(state_mrca,edge[,1:ncol(mat)+3])
  state_ans<-state_ans[order(c(nrow(mat)+1,phy$edge[,2])),]
  
  return(list(phy=phy,ancestral=state_ans))
}

greedy.partial.bit<-function(inmat,outmat=NULL,keep_cutoff=1,weight=NULL,noise=F,file_backed=NULL,verbose=F){
  phy<-greedy.wrap.bit(inmat,weight=weight,noise=noise,file_backed=file_backed,verbose=verbose)
  if(is.null(outmat)){
    return(phy$phy)
  }
  node_mat<-matrix(phy$ancestral[(1:nrow(inmat))*(-1),],,ncol(inmat))
  phy<-phy$phy
  outmat_retain<-apply(outmat,1,function(x){
    mutationNestedBit(x,node_mat[1,])
  })
  outmat_retain<-outmat_retain<keep_cutoff
  if(sum(outmat_retain)==0){
    return(list(phys=list(phy),category=rep(1,nrow(inmat)),isbroken=F))
  }else if(sum(outmat_retain)==1){
    outmat<-matrix(outmat[outmat_retain,],,ncol(inmat))
  }else{
    outmat<-outmat[outmat_retain,]
  }
  
  nodes_retain<-nodeKeepBit(edge=phy$edge,
                            Ntip=Ntip(phy),
                            node_mat=node_mat,
                            outmat=outmat,
                            keep_cutoff=keep_cutoff)
  nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
  broken_phys<-sub.multi.tree(phy,nodes_retain)
  
  category<-rep(0,nrow(inmat))
  i<-1
  while(i<=length(nodes_retain)){
    temp<-edge.subtree(phy,nodes_retain[i])
    category[temp[temp[,2]<=Ntip(phy),2]]<-i
    i<-i+1
  }
  
  names(category)<-rownames(inmat)
  return(list(phys=broken_phys,category=category,isbroken=T,origin_phy=phy,nodes_retain=nodes_retain))
}

split.mg.big<-function(binary,maxChar,weight=NULL,radius,p_shared=0.8,sqrt_nNeighbor=F,min_size=10,batch_size=1e4,file_backed=NULL,ncore=1,verbose=F){
  if(is.null(weight)){
    neighbor.fun<-function(i){
      radiusUnweightedParallel(
        i=i-1,mat=binary,maxChar=maxChar,
        radius=radius
      )
    }
  }else{
    neighbor.fun<-function(i){
      radiusWeightedParallel(
        i=i-1,mat=binary,maxChar=maxChar,
        weight=weight,radius=radius
      )
    }
  }
  if(!is.null(file_backed)){
    if(file_backed==""){
      neighbor_mat<-big.matrix(
        nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),
        type="integer",backingfile=""
      )
      backfile<-file.name(neighbor_mat)
    }else{
      fpath<-gsub("[^/]+$","",file_backed)
      fname<-gsub(fpath,"",file_backed)
      neighbor_mat<-big.matrix(
        nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),type="integer",
        backingfile=fname,backingpath=fpath,descriptorfile=paste0(fname,".desc")
      )
      backfile<-file_backed
    }

    if(verbose){
      cat("backup file:",backfile,fill=T)
    }
  }else{
    neighbor_mat<-big.matrix(nrow=ceiling(nrow(binary)/32),ncol=nrow(binary),type="integer",backingfile=NULL)
  }
  nNeighbor<-integer()
  for(i in 1:ceiling(nrow(binary)/batch_size)-1){
    offset<-i*batch_size
    row_process<-(offset+1):min(offset+batch_size,nrow(binary))
    if(ncore>1){
      cl<-makeForkCluster(ncore)
      neighbor_tmp<-parSapply(cl,row_process,neighbor.fun)
      stopCluster(cl)
      rm(cl)
    }else{
      neighbor_tmp<-sapply(row_process,neighbor.fun)
    }
    neighbor_mat[,row_process]<-neighbor_tmp[-1,]
    nNeighbor<-c(nNeighbor,neighbor_tmp[1,])
    rm(neighbor_tmp)
    gc()
    if(verbose){
      cat(as.character(Sys.time()),": batch",i+1,"finished",fill=T)
    }
  }
  rm(binary,weight)
  gc()
  temp<-tipMergeNeighborBig(neighbor_mat@address,nNeighbor=nNeighbor,
                            p=p_shared,sqrt_nNeighbor=sqrt_nNeighbor,verbose=verbose)
  temp<-regularize.split(temp,min_size=min_size)
  rm(neighbor_mat)
  gc()
  if(!is.null(file_backed)){
    if(file.exists(backfile)){
      unlink(backfile)
    }
    if(file.exists(paste0(backfile,".desc"))){
      unlink(paste0(backfile,".desc"))
    }
  }
  return(temp)
}

#' perform iterative imputation to fill the missing entries.
#'
#' @param mat The mutation matrix.
#' @param weight Weights of indels.
#' @param n_neighbor The number of neighbors used.
#' @param p_mutation The minimum percentage of shared mutations for neighbors.
#' @param n_iter The maximum number of iterations.
#' @param cutoff If the number of newly imputed entries is smaller than this value, the function will change the number of neighbors according to `n_neighbor`.
#' @param verobse Whether to print detailed information.
#' @param ncore The number of threads used.
#' 
#' @examples
#' \dontrun{
#' imputed<-iter.impute.packed(mat)
#'}
#'
#' @return an imputed matrix.

iter.impute.packed<-function(mat,weight=NULL,n_neighbor=10:2,p_mutation=0.5,n_iter=100,cutoff=max(floor(nrow(mat)/100),100),verbose=F,ncore=1){
  binary<-mat2bit(mat)
  nloci<-rev(binary$loci)[1]
  
  if (!is.null(weight)&class(weight)[1]=="matrix") {
    temp <- diff(binary$loci)
    weight <- unlist(lapply(which(temp != 0), function(i) {
      weight[1:temp[i], i]
    }))
  }
  
  id<-rownames(mat)
  maxChar<-ncol(mat)
  
  n_mask<-sum(mat==-1)
  total_imputed<-0
  k<-0
  neighbor_idx<-1
  initial<-T
  reinitial<-0
  
  if(is.null(weight)){
    fun1<-function(i){
      findNeighborUnweightedParallel(i=i-1,mat=binary,maxChar=maxChar,
                                     n=n_neighbor[neighbor_idx],p=p_mutation)
    }
    fun2<-function(i){
      findNeighborRangeUnweighted(i=i-1,mat=binary,maxChar=maxChar,
                                  range=range_found,offset=0,
                                  n=n_neighbor[neighbor_idx],p=p_mutation)
    }
  }else{
    fun1<-function(i){
      findNeighborWeightedParallel(i=i-1,mat=binary,maxChar=maxChar,
                                   weight=weight,n=n_neighbor[neighbor_idx],p=p_mutation)
    }
    fun2<-function(i){
      findNeighborRangeWeighted(i=i-1,mat=binary,maxChar=maxChar,
                                range=range_found,offset=0,
                                weight=weight,n=n_neighbor[neighbor_idx],p=p_mutation)
    }
  }
  impute.fun1<-function(i){
    imputeMat(i-1,mat,range_found,0)
  }
  impute.fun2<-function(i){
    imputeMat(i-1,mat,neighbor_found,0)
  }  
  {
    if(verbose){
      cat(as.character(Sys.time()),": iteration ",k+1,", No. neighbor ",n_neighbor[neighbor_idx],fill=T)
    }
    binary<-mat2bit(mat)$binary
    if(ncore>1){
      cl<-makeForkCluster(ncore)
      neighbor_found<-t(parSapply(cl,1:nrow(binary),fun1))
      mat<-t(parSapply(cl,1:nrow(mat),impute.fun2))
      stopCluster(cl)
    }else{
      neighbor_found<-t(sapply(1:nrow(binary),fun1))
      mat<-t(sapply(1:nrow(mat),impute.fun2))
    }
    k<-k+1
    temp<-n_mask-sum(mat==-1)
    total_imputed<-total_imputed+temp
    if(verbose){
      cat("Impute ",temp," entries.",fill=T)
    }
    n_mask<-n_mask-temp
  }
  gc()
  while(k<n_iter&neighbor_idx<=length(n_neighbor)){
    if(verbose){
      cat(as.character(Sys.time()),": iteration ",k+1,", No. neighbor ",n_neighbor[neighbor_idx],fill=T)
    }
    binary<-mat2bit(mat)$binary
    if(initial){
      if(verbose){
        cat("initializing",fill=T)
      }
      if(ncore>1){
        cl<-makeForkCluster(ncore)
        range_found<-t(parSapply(cl,1:nrow(binary),fun1))
        mat<-t(parSapply(cl,1:nrow(mat),impute.fun1))
        stopCluster(cl)
      }else{
        range_found<-t(sapply(1:nrow(binary),fun1))
        mat<-t(sapply(1:nrow(mat),impute.fun1))
      }
      initial<-F
    }else{
      if(ncore>1){
        cl<-makeForkCluster(ncore)
        neighbor_found<-t(parSapply(cl,1:nrow(mat),fun2))
        mat<-t(parSapply(cl,1:nrow(mat),impute.fun2))
        stopCluster(cl)
      }else{
        neighbor_found<-t(sapply(1:nrow(mat),fun2))
        mat<-t(sapply(1:nrow(mat),impute.fun2))
      }
    }
    gc()
    k<-k+1
    temp<-n_mask-sum(mat==-1)
    total_imputed<-total_imputed+temp
    if(verbose){
      cat("Impute ",temp," entries.",fill=T)
    }
    n_mask<-n_mask-temp
    if(temp<cutoff){
      neighbor_idx<-neighbor_idx+1
      reinitial<-reinitial+1
      if(reinitial>5){
        initial<-T
        reinitial<-0
      }
    }
  }
  rownames(mat)<-id
  if(verbose){
    cat("Impute ",k," times with ",total_imputed," entries in total.",sep="",fill=T)
  }
  return(mat)
}

write.mat2bin<-function(mat,file){
  output<-rev(dim(mat))
  output<-c(output,mat)
  writeBin(object=output,con=file)
}


#' read packed matrix from bin file
#'
#' @param file The file name.
#' @param chunk_size The chunk size used in reading.
#' 
#' @examples
#' \dontrun{
#' mat<-read.bin2mat(file="test.bin")
#'}
#'
#' @return an bit coded matrix containing the mutation profiles.

read.bin2mat<-function(file,chunk_size=32768L){
  con<-file(file,"rb")
  d <- readBin(con = con, what = "integer", n = 2)
  if(log2(d[1])+log2(d[2])<31){
    mat <- readBin(con = con, what = "integer", n = d[1] * 
                     d[2] + 2)
    mat <- matrix(mat, d[2], d[1])
  }else{
    mat<-matrix(integer(),d[2],d[1])
    i<-0L
    while(i<d[1]){
      chunk_data<-readBin(con, what = "integer", n = d[2]*chunk_size)
      mat[,(1+i):min((chunk_size+i),d[1])]<-chunk_data
      i<-i+chunk_size
    }
  }
  close(con)
  return(mat)
}

bs.collapse<-function(tre,bs,threshold=50){
  len<-bsCollapse(tre$edge,tre$edge.length,bs,threshold)
  tre$edge.length<-len
  tre<-di2multi(tre)
  return(tre)
}

phy.join<-function(
    inmat,bit_coded=F,weight=NULL,ancestralType=0,
    tipSize=rep(1L,ncol(inmat)),minSize=0,minRatio=2
){
  phy<-list()
  ntip<-ncol(inmat)
  tipLabel<-colnames(inmat)
  if(bit_coded){
    tmp<-greedyJoinEdgeBitReverse(inmat,weight=weight,verbose=F,ancestralType=ancestralType)
    rm(inmat)
    nodeStateIdx<-c(tmp$edge[2,]+1L,ntip*2-1)
    edgeTmp<-tmp$edge+1L
    edge<-reOrderEdge2(edgeTmp,ntip*2-1)
    uniqNode<-unique(as.vector(edge))
    phy$edge<-t(matrix(as.integer(factor(edge,levels=c(1:ntip,uniqNode[uniqNode>ntip]))),2,))
    transition<-rep(0,ntip*2-1)
    transition[edge[2,]]<-phy$edge[,2]
    transition[ntip*2-1]<-ntip+1
    nodeState<-tmp$nodeState[,order(transition[nodeStateIdx])]
    nodeBistate<-tmp$nodeBistate[,order(transition[nodeStateIdx])]
    phy$Nnode<-max(phy$edge)-ntip
    phy$tip.label<-tipLabel
    phy$edge.length <- rep(0, nrow(phy$edge))
    class(phy)<-"phylo"
    attr(phy,"order")<-"cladewise"
    return(list(phy=phy,node_mat=nodeState[,1:(ntip-1)+ntip],nodeBistate=nodeBistate))
  }else{
    edge<-greedyJoinEdge2(inmat,weight=weight,tipSize=tipSize,minSize=minSize,minRatio=minRatio,verbose=F)
    state_mrca<-edge[1:nrow(inmat)+3,ncol(edge)]
    edge<-edge[,(-1)*ncol(edge)]
    edge[1:2,]<-edge[1:2,]+1
    edge<-reOrderEdge2(edge,ncol(inmat)*2-1)
    phy$edge.length<-edge[3,]
    phy$edge<-t(edge[1:2,])
    temp<-unique(as.vector(phy$edge))
    temp<-temp[temp>ncol(inmat)]
    phy$edge<-matrix(as.integer(factor(phy$edge,levels=c(1:ncol(inmat),temp))),,2)
    phy$Nnode<-max(phy$edge)-ncol(inmat)
    phy$tip.label<-colnames(inmat)
    class(phy)<-"phylo"
    attr(phy,"order")<-"cladewise"
    
    node_mat<-cbind(state_mrca,edge[1:nrow(inmat)+3,phy$edge[,2]>ncol(inmat)])[,order(c(ncol(inmat)+1,phy$edge[phy$edge[,2]>ncol(inmat),2]))]
    return(list(phy=phy,node_mat=node_mat))
  }
}

edge2phy<-function(edge,tipLabel){
  phy<-list()
  phy$edge<-t(edge[1:2,]+1)
  phy$edge.length<-rep(0,nrow(phy$edge))
  phy$Nnode<-length(tipLabel)-1
  phy$tip.label<-tipLabel
  class(phy)<-"phylo"
  attr(phy,"order")<-"cladewise"
  return(phy)
}

phy.break<-function(phy,nodes){
  tmp<-edgeSplit(phy$edge,nodes)+1
  category<-tmp[match(1:Ntip(phy),phy$edge[,2])]
  names(category)<-phy$tip.label
  phys<-lapply(1:length(nodes),function(i){
    subEdge<-phy$edge[tmp==i,]
    output<-list()
    output$tip.label<-phy$tip.label[sort(subEdge[subEdge[,2]<=Ntip(phy),2])]
    output$edge<-matrix(as.integer(factor(subEdge)),,2)
    output$edge.length<-phy$edge.length[tmp==i]
    output$Nnode<-max(output$edge)-length(output$tip.label)
    class(output)<-"phylo"
    attr(output,"order")<-"cladewise"
    return(output)
  })
  return(list(phys=phys,category=category))
}

member.combine<-function(subMember){
  offset<-diffinv(sapply(subMember,function(x) max(x)+1))
  unlist(lapply(1:length(subMember),function(i){
    tmp<-subMember[[i]]
    tmp[tmp!=-1]+offset[i]
  }))
}

#' The top solver that splits the DNA sequence data into small clusters
#'
#' @param mat The mutation matrix.
#' @param splt The result from previous run.
#' @param split_idx The index of cluster to further split.
#' @param radius The minimum number of mutations between neighbors.
#' @param p_shared The threshold for the percentage of shared neighbors.
#' @param min_size The minimum size of clusters. Clusters smaller than this value will be collapsed.
#' @param file_backed The name of backup file. If set to `NULL`, do not create backup file. If set to `""`, create an anonymous backup file.
#' @param ncore The number of threads used.
#' @param cuda_path The path of CUDA script. If set to `NULL`, do not use CUDA.
#' @param verbose Whether print detailed information.
#' 
#' 
#' @examples
#' \dontrun{
#' splt<-split.data.DNA(mat,radius=1000,p_shared=0.8)
#'}
#'
#' @return Return an integer matrix containing the cluster information.

split.data.DNA<-function(mat,splt=NULL,splt_idx=1,
                         radius=100,p_shared=0.8,min_size=10,
                         file_backed=NULL,ncore=1,cuda_path=NULL,verbose=F){
  batch_size=1e4
  if(is.null(weight)){
    weight<-rep(1,nrow(mat)*4)
  }
  
  if(is.null(splt)){
    layer<-1
    k<-1
    splt<-cbind(rep(1,ncol(mat)),rep(0,ncol(mat)))
  }else{
    layer<-which(apply(splt,2,function(x){
      splt_idx%in%x
    }))
    k<-max(splt)
    if(layer==ncol(splt)){
      splt<-cbind(splt,rep(0,ncol(mat)))
    }
  }
  
  idx<-which(splt[,layer]==splt_idx)
  mat<-mat[,idx]
  
  if(is.null(cuda_path)){
    neighbor.fun<-function(i){
      radiusWeightedDNAParallel(
        i=i-1,mat=mat,
        weight=weight,radius=radius
      )
    }
    if(!is.null(file_backed)){
      if(file_backed==""){
        neighbor_mat<-big.matrix(
          nrow=ceiling(ncol(mat)/32),ncol=ncol(mat),
          type="integer",backingfile=""
        )
        backfile<-file.name(neighbor_mat)
      }else{
        fpath<-gsub("[^/]+$","",file_backed)
        fname<-paste0(gsub(fpath,"",file_backed),".bk")
        neighbor_mat<-big.matrix(
          nrow=ceiling(ncol(mat)/32),ncol=ncol(mat),type="integer",
          backingfile=fname,backingpath=fpath,descriptorfile=paste0(fname,".desc")
        )
        backfile<-paste0(file_backed,".bk")
      }
      if(verbose){
        cat("backup file:",backfile,fill=T)
      }
    }else{
      neighbor_mat<-big.matrix(
        nrow=ceiling(ncol(mat)/32),ncol=ncol(mat),
        type="integer",backingfile=NULL
      )
    }
    nNeighbor<-integer()
    for(i in 1:ceiling(ncol(mat)/batch_size)-1){
      offset<-i*batch_size
      row_process<-(offset+1):min(offset+batch_size,ncol(mat))
      if(ncore>1){
        cl<-makeForkCluster(ncore)
        neighbor_tmp<-parSapply(cl,row_process,neighbor.fun)
        stopCluster(cl)
        rm(cl)
      }else{
        neighbor_tmp<-sapply(row_process,neighbor.fun)
      }
      neighbor_mat[,row_process]<-neighbor_tmp[-1,]
      nNeighbor<-c(nNeighbor,neighbor_tmp[1,])
      rm(neighbor_tmp)
      gc()
      if(verbose){
        cat(as.character(Sys.time()),": batch",i+1,"finished",fill=T)
      }
    }
    rm(mat,weight)
    gc()
    ncoreMerge<-min(ncore,ceiling(length(idx)/2e4))
    if(verbose){
      cat(as.character(Sys.time()),": start tipMerge with ",ncoreMerge," threads.",sep="",fill=T)
    }
    if(length(idx)<=batch_size||ncoreMerge==1||ncoreMerge==2){
      temp<-tipMergeNeighborBig(neighbor_mat@address,nNeighbor=nNeighbor,
                                p=p_shared,sqrt_nNeighbor=F,verbose=verbose)
    }else{
      nbatch<-ncoreMerge
      bsize<-ceiling(ncol(neighbor_mat)/nbatch)
      desc<-describe(neighbor_mat)
      cl<-makeForkCluster(ncoreMerge)
      mb<-parLapply(cl,1:nbatch-1,function(i){
        bgmat<-attach.big.matrix(desc)
        tipMergeSub(
          bgmat@address,
          i*bsize,bsize,
          nNeighbor,
          p=p_shared
        )
      })
      stopCluster(cl)
      if(verbose){
        cat("=")
      }
      mb<-member.combine(mb)
      tmp<-matrix(,nbatch,nbatch)
      colIdx<-col(tmp)[which(lower.tri(tmp))]-1
      rowIdx<-row(tmp)[which(lower.tri(tmp))]-1
      trans<-0:max(mb)
      offset<-0
      ntotal<-(nbatch*(nbatch-1)/2)
      while(offset<ntotal){
        processIdx<-(offset+1):min(offset+nbatch,ntotal)
        cl<-makeForkCluster(ncoreMerge)
        transTmp<-parSapply(cl,processIdx,function(i){
          bgmat<-attach.big.matrix(desc)
          offset1<-(colIdx[i])*bsize
          offset2<-(rowIdx[i])*bsize
          res<-tipMergeCross(
            bgmat@address,
            offset1,offset2,bsize,
            nNeighbor,mb,
            trans,p=p_shared
          )
          return(res)
        })
        stopCluster(cl)
        while(ncol(transTmp)>10){
          nRes<-ncol(transTmp)
          cl<-makeForkCluster(floor(nRes/2))
          if(nRes%%2==0){
            res<-parApply(cl,matrix(transTmp,,nRes/2),2,function(x){
              x<-matrix(x,,2)
              return(combineTransition(x[,1],x[,2]))
            })
            transTmp<-res
          }else{
            res<-parApply(cl,matrix(transTmp[,1:(nRes-1)],,(nRes-1)/2),2,function(x){
              x<-matrix(x,,2)
              return(combineTransition(x[,1],x[,2]))
            })
            transTmp<-cbind(res,transTmp[,nRes])
          }
          stopCluster(cl)
        }
        trans<-transTmp[,1]
        for(i in 2:ncol(transTmp)){
          trans<-combineTransition(trans,transTmp[,i])
        }
        offset<-offset+nbatch
        if(verbose){
          cat("=")        
        }
      }
      if(verbose){
        cat("",fill=T)
      }
      temp<-sapply(mb,function(i){
        transitionRoot(trans,i)
      })
    }
    rm(neighbor_mat)
    gc()
    if(verbose){
      cat(as.character(Sys.time()),": tipMerge finished.",sep="",fill=T)
    }
    if(!is.null(file_backed)){
      if(file.exists(backfile)){
        unlink(backfile)
      }
      if(file.exists(paste0(backfile,".desc"))){
        unlink(paste0(backfile,".desc"))
      }
    }
  }else{
    if(is.null(file_backed)||file_backed==""){
      cat("file_backed must be specified when using cuda.",fill=T)
      return(NULL)
    }
    write.mat2bin(mat,file=paste0(file_backed,".tbmat"))
    writeBin(object=as.vector(weight),con=paste0(file_backed,".w"))
    rm(mat)
    cat(as.character(Sys.time()),": CUDA start.",sep="",fill=T)
    system(paste0(cuda_path,"radiusNeighborDNA ",file_backed," ",radius))
    system(paste0(cuda_path,"tipMerge ",file_backed," ",p_shared))
    cat(as.character(Sys.time()),": CUDA stop.",sep="",fill=T)
    unlink(paste0(file_backed,".tbmat"))
    unlink(paste0(file_backed,".w"))
    unlink(paste0(file_backed,".bk"))
    unlink(paste0(file_backed,".nNeighbor"))
    unlink(paste0(file_backed,".idx"))
    temp<-read.table(file=paste0(file_backed,".member"))[,1]
    unlink(paste0(file_backed,".member"))
  }
  temp<-regularize.split(temp,min_size=min_size)
  temp[temp!=0]<-temp[temp!=0]+k
  splt[idx,layer+1]<-temp
  return(splt)
}

write.mat2fa<-function(mat,file,randomize_ambiguous=F){
  seq_name<-colnames(mat)
  decoded<-apply(mat,2,decodeATGC8,randomize_ambiguous=randomize_ambiguous)
  output<-paste0(">",colnames(mat),"\n",decoded)
  cat(output,file=file,sep="\n")
}

get.break.node<-function(rootVote,cutoff=200,pct_threshold=0.05){
  retain<-which(rootVote[,2]>=max(rootVote[,2])-cutoff)
  if(length(retain)==0){
    tmp<-table(rootVote[,1])/nrow(rootVote)
  }else{
    tmp<-table(rootVote[retain,1])/length(retain)
  }
  root_node<-as.integer(names(tmp)[which.max(tmp)])
  break_node<-setdiff(as.integer(names(tmp)[tmp>=pct_threshold]),root_node)
  return(list(root_node=root_node,break_node=break_node))
}

break.pp<-function(phy,root_node,break_node){
  phy$node.label<-1:Nnode(phy)+Ntip(phy)
  break_edge<-matrix(phy$edge[phy$edge[,2]%in%break_node,],,2)
  if(root_node>Ntip(phy)){
    phy2<-root(phy,node=root_node,resolve.root=T)
  }else{
    phy2<-root(phy,outgroup=root_node,resolve.root=T)
  }
  transition<-rep(0,Nnode(phy2))
  transitionIdx<-unique(phy2$edge[phy2$edge>Ntip(phy2)])-Ntip(phy2)
  transition[transitionIdx]<-1:Nnode(phy2)+Ntip(phy2) ##old idx to new idx
  phy2$node.label[transition[1:Nnode(phy2)]-Ntip(phy2)]<-phy2$node.label
  phy2$edge[phy2$edge>Ntip(phy2)]<-transition[phy2$edge[phy2$edge>Ntip(phy2)]-Ntip(phy2)]
  if(length(break_node)==0){
    return(list(phy=phy2,break_edge=matrix(,0,2),break_node=integer()))
  }
  
  transition<-1:(Nnode(phy)+Ntip(phy))
  transition[as.integer(phy2$node.label[-1])]<-2:Nnode(phy2)+Ntip(phy2)
  break_edge2<-t(apply(break_edge,1,function(x){
    transition[x]
  }))
  phy2$node.label<-1:Nnode(phy2)+Ntip(phy2)
  break_node2<-apply(break_edge2,1,function(x){
    if(x[2]<=Ntip(phy2)){
      return(x[2])
    }else{
      return(max(x))
    }
  })
  return(list(phy=phy2,break_edge=break_edge2,break_node=break_node2))
}

gen.mat.big<-function(mat,tipname,ansMat=NULL,nodename=NULL,splt,return.mat=T){
  mg_idx<-setdiff(unique(splt),0)
  if(length(mg_idx)==0){
    submat_idx<-1:ncol(mat)
    if(return.mat){
      output<-mat
      colnames(output)<-tipname
      return(list(submat=output,idx=submat_idx,splt=rep(0,ncol(mat))))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=rep(0,ncol(mat))))
    }
  }else{
    if(return.mat){
      output<-cbind(
        output<-mat[,splt==0],
        ansMat[,match(paste0("__M",mg_idx),nodename)]
      )
      colnames(output)<-c(tipname[splt==0],paste0("__M",mg_idx))
    }
    offset<-sum(splt==0)
    submat_idx<-sapply(splt,function(x){
      if(x==0) return(0)
      else return(match(x,mg_idx)+offset)
    })
    submat_idx[submat_idx==0]<-1:offset
    if(return.mat){
      return(list(submat=output,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }else{
      return(list(submat=NULL,idx=submat_idx,splt=c(rep(0,sum(splt==0)),mg_idx)))
    }
  }
}

#' read packed matrix from bin file into a big matrix
#'
#' @param file The file name.
#' @param chunk_size The chunk size used in reading.
#' @param backingfile The file name of backing file. If NULL, create a big matrix in memory.
#' @param backingpath The path to the backing file.
#' @param descriptorfile The file name of descriptor file.
#' 
#' 
#' @examples
#' \dontrun{
#' mat<-read.bin2bigmat(file="test.bin")
#'}
#'
#' @return a big matrix.

read.bin2bigmat<-function(file,chunk_size=32768L,
                          backingfile=NULL,backingpath=NULL,descriptorfile=NULL){
  con<-file(file,"rb")
  d <- readBin(con = con, what = "integer", n = 2)
  mat<-big.matrix(
    nrow=d[2],ncol=d[1],type="integer",
    backingfile=backingfile,backingpath=backingpath,descriptorfile=descriptorfile
  )
  if(log2(d[1])+log2(d[2])<31){
    mat[1:(d[1]*d[2])] <- readBin(con = con, what = "integer", n = d[1] * 
                                    d[2] + 2)
  }else{
    i<-0L
    while(i<d[1]){
      chunk_data<-readBin(con, what = "integer", n = d[2]*chunk_size)
      mat[,(1+i):min((chunk_size+i),d[1])]<-chunk_data
      i<-i+chunk_size
    }
  }
  close(con)
  return(mat)
}

#' reconstruct cluster-level phylogenies
#'
#' @param j The index of cluster.
#' @param mat_desc The descriptor of matrix.
#' @param tipname The names of tips.
#' @param matfile The file name of backing matrix in hard disk.
#' @param ansMat_desc The descriptor of matrix of pseudo-OTUs.
#' @param nodename The name of pseudo-OTUs.
#' @param splt Split of data from `split.data.DNA`.
#' @param layer The layer of clusters.
#' @param prefix The prefix of temporary file.
#' @param cuda_path The path of CUDA script. If set to `NULL`, do not use CUDA.
#' @param sample_size The percentage of outgroups sampled.
#' @param pct_threshold The frequency threshold used in breaking edges.
#' @param nbp_threshold The threshold used to filter out distal related outgroups.
#' @param vft_thread The number of threads used in VeryFastTree.
#' 
#' @examples
#' \dontrun{
#' res<-partial.join.DNA(
#'        mat_desc=mat_desc,
#'        tipname=tipname,
#'        matfile="test.bin",
#'        ansMat_desc=NULL,
#'        nodename=NULL,
#'        splt=splt,
#'        layer=2,
#'        prefix="test",
#'        cuda_path="~/CasPhy-cuda/build/"
#'      )
#'}
#'
#' @return a list containing the result.

partial.join.DNA<-function(j,mat_desc,tipname,matfile,ansMat_desc,nodename,splt,layer,prefix,cuda_path,sample_size=0.1,
                       pct_threshold=0.001,nbp_threshold=200,
                       vft_thread=8){
  mat<-attach.big.matrix(mat_desc)
  if(!is.null(ansMat_desc)){
    ansMat<-attach.big.matrix(ansMat_desc)
  }else{
    ansMat<-NULL
  }
  inIdx<-which(splt[,layer]==j)
  outIdx<-which(splt[,layer]!=j)
  inmat<-gen.mat.big(mat[,inIdx],tipname[inIdx],ansMat,nodename,splt[inIdx,layer+1])
  tip2pseudo<-inmat$idx
  inmat<-inmat$submat
  tipLabel<-colnames(inmat)
  
  ##run VeryFastTree
  write.mat2fa(inmat,file=paste0(prefix,"__M",j,".fasta"),randomize_ambiguous=F)
  system(paste0("VeryFastTree -nosupport -fastest -nt ",prefix,"__M",j,".fasta -threads ",vft_thread," 1>",prefix,"__M",j,".nwk 2>",prefix,"__M",j,".vftlog"))
  unlink(paste0(prefix,"__M",j,".fasta",sep=""))
  phy<-read.tree(paste0(prefix,"__M",j,".nwk"))
  unlink(paste0(prefix,"__M",j,".nwk"))
  phy<-order.tip(phy,tipLabel)
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-cbind(inmat,node_mat)
  t2<-Sys.time()
  
  if(length(outIdx)==0){
    return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
  }
  nSample<-max(16384,floor(length(outIdx)*sample_size))
  if(length(outIdx)>nSample){
    outIdx<-sort(sample(outIdx,nSample))
  }
  
  ##root phy
  t1<-Sys.time()
  if(Ntip(phy)<=100*16384/nSample){
    #    rootVote<-findRoot(phy$edge,node_mat,mat,outIdx-1L)
    rootVote<-findRootAVX2Big(phy$edge,node_mat,mat@address,outIdx-1L)
  }else{
    write.mat2bin(node_mat,file=paste0(prefix,"__M",j,".nodeState"))
    write.mat2bin(matrix(outIdx-1L,,1),file=paste0(prefix,"__M",j,".outmatIdx"))
    write.mat2bin(t(phy$edge),file=paste0(prefix,"__M",j,".edge"))
    system(paste0(cuda_path,"voteRoot2 ",matfile," ",prefix,"__M",j))
    rootVote<-read.bin2mat(paste0(prefix,"__M",j,".rootVote"))
    rootVote[,1]<-rootVote[,1]+1L
    unlink(paste0(prefix,"__M",j,".nodeState"))
    unlink(paste0(prefix,"__M",j,".outmatIdx"))
    unlink(paste0(prefix,"__M",j,".edge"))
    unlink(paste0(prefix,"__M",j,".rootVote"))
  }
  breakRes<-get.break.node(rootVote,cutoff=nbp_threshold,pct_threshold=pct_threshold)
  root_node<-breakRes$root_node
  break_node<-breakRes$break_node
  rootedRes<-break.pp(phy,root_node,break_node)
  phy<-rootedRes$phy
  break_node<-rootedRes$break_node
  nodes_retain<-breakUnrootedPhy(phy$edge,Ntip(phy),Nnode(phy),breakNode=break_node)
  nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
  if(length(nodes_retain)==0){
    return(list(phys=list(),category=rep(0,length(inIdx)),ansState=NULL,isbroken=T,inIdx=inIdx))
  }
  t2<-Sys.time()
  
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-node_mat[,nodes_retain-Ntip(phy)]
  tmp<-phy.break(phy,nodes_retain)
  broken_phys<-tmp$phys
  category<-tmp$category
  return(list(phys=broken_phys,category=category[tip2pseudo],ansState=node_mat,isbroken=T,inIdx=inIdx))
}

#' summarize the results from `partial.join.DNA`
#'
#' @param partial_res_collected The results from `partial.join.DNA`.
#' @param clustIdx The indices of clusters.
#' @param splt The splt object used in phylogenetic reconstruction.
#' @param layer The layer of clusters.
#' 
#' @examples
#' \dontrun{
#' mat<-read.bin2bigmat(file="test.bin")
#'}
#'
#' @return summarized results.

layer.summarize<-function(partial_res_collected,clustIdx,splt,layer){
  max_mg_idx<-max(splt)
  nPhys<-sapply(partial_res_collected,function(x){
    if(x$isbroken){
      return(length(x$phys))
    }else{
      return(0)
    }
  })
  offset<-diffinv(nPhys)+max_mg_idx
  tmp<-unlist(lapply(1:length(clustIdx),function(i){
    if(partial_res_collected[[i]]$isbroken){
      if(nPhys[i]==0){
        return(integer())
      }else{
        return(1:nPhys[i]+offset[i])
      }
    }else{
      return(clustIdx[i])
    }
  }))
  newPhys<-unlist(lapply(partial_res_collected,function(x) x$phys),recursive=F)
  if(length(tmp)>0){
    names(newPhys)<-paste0("__M",tmp)
  }
  # phys<-c(phys,newPhys)
  
  newAnsMat<-lapply(partial_res_collected,function(x) x$ansState)
  nbin<-NULL
  for(obj in newAnsMat){
    if(!is.null(obj)){
      if(is.matrix(obj)){
        nbin<-nrow(obj)
      }else{
        nbin<-length(obj)
      }
    }
  }
  if(!is.null(nbin)){
    newAnsMat<-matrix(unlist(newAnsMat),nbin,)
  }else{
    newAnsMat<-NULL
  }
  if(length(tmp)>0&!is.null(newAnsMat)){
    colnames(newAnsMat)<-paste0("__M",tmp)      
  }
  # ansMat<-cbind(ansMat,newAnsMat)
  
  replaceIdx<-unlist(lapply(partial_res_collected,function(x){
    if(x$isbroken){
      return(x$inIdx)
    }else{
      return(integer())
    }
  }))
  new_group<-unlist(lapply(1:length(clustIdx),function(i){
    if(partial_res_collected[[i]]$isbroken){
      tmp<-partial_res_collected[[i]]$category
      tmp[tmp!=0]<-tmp[tmp!=0]+offset[i]
      return(tmp)
    }else{
      return(integer())
    }
  }))
  temp<-splt[replaceIdx,layer+1]
  temp[new_group>0]<-new_group[new_group>0]
  splt[replaceIdx,layer]<-temp
  return(list(splt=splt,phys=newPhys,ansMat=newAnsMat))
}

