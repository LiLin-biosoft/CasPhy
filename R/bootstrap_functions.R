phy2clade<-function(phy,prefix){
  ntip<-Ntip(phy)
  node_size<-nodeDesBig(phy$edge,Nnode(phy),ntip,file_path=paste0(prefix,".bk"),header=F)
  
  fpath<-gsub("[^/]+$","",prefix)
  if(fpath==""){
    fpath="./"
    fname=prefix
  }else{
    fname<-gsub(fpath,"",prefix)
  }
  
  ncols<-Nnode(phy)-1L
  nrows<-as.integer(ceiling(ntip/32))
  desc<-new("big.matrix.descriptor")
  desc@description<-list(
    sharedType="FileBacked",
    filename=paste0(fname,".bk"),
    dirname=fpath,
    totalRows=nrows,
    totalCols=ncols,
    rowOffset=c(0,nrows),
    colOffset=c(0,ncols),
    nrow=nrows+0,
    ncol=ncols+0,
    rowNames=NULL,
    colNames=NULL,
    type="integer",
    separated=F
  )
  return(list(node_size=c(ntip,rev(node_size)),desc=desc,filename=paste0(prefix,".bk")))
}

clade.support<-function(phys,prefix,threshold=0.05,ncore=10){
  compare.fun<-function(x){
    bitClade<-attach.big.matrix(desc)
    resCmp<-approximateCompare2_AVX2(bitClade@address,ntip,phyRef$edge,x$edge,Nnode(x),threshold=threshold+0.05)
    rm(bitClade)
    tmp<-gc()
    return(resCmp)
  }
  ntip<-Ntip(phys[[1]])
  output<-list()
  for(i in 1:ceiling(length(phys)/2)){
    cat("Run",i,fill=T)
    phyRef<-phys[[i]]
    clade_mat<-phy2clade(phyRef,prefix)
    desc<-clade_mat$desc
    
    if(ncore>1){
      cl<-makeForkCluster(min(ncore,length(phys)-i))
      res<-parLapply(cl,phys[(i+1):length(phys)],compare.fun)
      stopCluster(cl)
    }else{
      res<-lapply(phys[(i+1):length(phys)],compare.fun)
    }
    
    unlink(clade_mat$filename)
    nodesQry<-sapply(res,function(x){
      x$ref2qry[,1]
    })
    nodesQry<-cbind(1:Nnode(phyRef)+ntip,nodesQry)
    xorDis<-sapply(res,function(x){
      x$ref2qry[,2]
    })
    xorDis<-cbind(0,xorDis)
    if(i>1){
      tmp<-lapply(output[1:(i-1)],function(x){
        qry2ref(x[[1]][,i],x[[2]][,i],ntip,Nnode(phys[[i]]))
      })
      nodesQry<-cbind(sapply(tmp,function(x) x[,1]),nodesQry)
      xorDis<-cbind(sapply(tmp,function(x) x[,2]),xorDis)
    }
    
    output[[i]]<-list(nodes=nodesQry,dis=xorDis,lib=clade_mat$node_size)
  }
  return(output)
}

mark.overlap<-function(nodeSets,nodeSize,ncore=1){
  if(length(nodeSets)>1){
    setSize<-sapply(nodeSets,function(x){
      sum(nodeSize[x])
    })
    odr<-order(setSize)
    nodeSets<-nodeSets[order(setSize)]
    ncoreRun<-min(ncore,length(nodeSets)-1)
    if(ncoreRun>1){
      cl<-makeForkCluster(ncoreRun)
      shufIdx<-sample(2:length(nodeSets)-1,length(nodeSets)-1)
      overlapped<-parSapply(cl,shufIdx,function(j){
        tmp<-sapply((j+1):length(nodeSets),function(k){
          length(intersect(nodeSets[[j]],nodeSets[[k]]))>0
        })
        return(sum(tmp)>0)
      })
      stopCluster(cl)
      overlapped<-overlapped[order(shufIdx)]
    }else{
      overlapped<-sapply(2:length(nodeSets)-1,function(j){
        tmp<-sapply((j+1):length(nodeSets),function(k){
          length(intersect(nodeSets[[j]],nodeSets[[k]]))>0
        })
        return(sum(tmp)>0)
      })
    }
    overlapped<-c(overlapped,F)
    return(overlapped[match(1:length(nodeSets),odr)])
  }else if(length(nodeSets)==1){
    return(F)
  }else{
    return(logical())
  }
}

node.relabel<-function(edge,ntip){
  uniqNode<-unique(as.vector(edge))
  edge<-t(matrix(as.integer(factor(edge,levels=c(1:ntip,uniqNode[uniqNode>ntip]))),2,))
  return(edge)
}

rebuild.nodes<-function(targetEdge,targetNodeSize,ntip,nodeClean,prefix,CUDA_path,threshold=0.05,ncore=10,check=T){
  maxIdx<-length(targetNodeSize)+1L
  while(length(nodeClean)>0){
    cat("No. total sets:",length(nodeClean),fill=T)
    ncoreRun<-min(ncore,length(nodeClean))
    if(ncoreRun>1){
      cl<-makeForkCluster(ncoreRun)
      tmpJoin<-parLapply(cl,nodeClean,function(x){
        buildNode(targetEdge,targetNodeSize,x,ntip,threshold=threshold)
      })
      stopCluster(cl)
    }else{
      tmpJoin<-lapply(nodeClean,function(x){
        buildNode(targetEdge,targetNodeSize,x,ntip,threshold=threshold)
      })
    }
    nJoin<-sapply(tmpJoin,length)
    cat("No. failed join:",sum(nJoin==0),fill=T)
    cat("No. redundant join:",sum(nJoin==1),fill=T)
    retain<-which(sapply(tmpJoin,length)>1)
    nodeClean<-nodeClean[retain]
    tmpJoin<-tmpJoin[retain]
    overlapped<-mark.overlap(tmpJoin,targetNodeSize,ncore=ncore)
    tmpJoin<-tmpJoin[!overlapped]
    nodeClean<-nodeClean[overlapped]
    cat("No. processing sets:",length(tmpJoin),fill=T)
    cat("No. remaining sets:",length(nodeClean),fill=T)
    cat("",fill=T)
    
    i<-1
    while(i<=length(tmpJoin)){
      edgeIdx<-which(targetEdge[,2]%in%tmpJoin[[i]])
      ansIdx<-targetEdge[edgeIdx[1],1]
      targetEdge[edgeIdx,1]<-maxIdx
      targetEdge<-rbind(targetEdge,c(ansIdx,maxIdx))
      targetNodeSize<-c(targetNodeSize,sum(targetNodeSize[tmpJoin[[i]]]))
      maxIdx=maxIdx+1L
      i<-i+1
    }
  }
  
  write.mat2bin(t(targetEdge)-1L,file=paste0(prefix,".edge"))
  cat("Run cmd:",paste0("edge2cladewise ",prefix," ",ntip," ",max(targetEdge)),fill=T)
  system(paste0(CUDA_path,"/edge2cladewise ",prefix," ",ntip," ",max(targetEdge)))
  tmpEdgeCUDA<-read.bin2mat(file=paste0(prefix,".edge2"))
  tmpEdgeCUDA<-t(tmpEdgeCUDA+1L)
  unlink(paste0(prefix,".edge2"))
  
  if(check){
    tmpEdge<-reOrderEdge2(t(targetEdge),ntip+1)
    tmpEdge<-node.relabel(tmpEdge,ntip)
    cat("No. difference between Rcpp and CUDA:",sum(tmpEdgeCUDA!=tmpEdge),fill=T)
  }
  return(tmpEdgeCUDA)
}

phy.wrap<-function(edge,tipLabel){
  ntip<-length(tipLabel)
  phy<-list()
  phy$tip.label<-tipLabel
  phy$edge<-edge
  phy$edge.length<-NULL
  phy$Nnode<-max(phy$edge)-ntip
  class(phy)<-"phylo"
  attr(phy,"order")<-"cladewise"
  return(phy)
}

cal.abs<-function(phy,bstrees,prefix,threshold=0.05,ncore=10){
  
  ntip<-Ntip(phy)
  clade_mat<-phy2clade(phy,prefix)
  desc<-clade_mat$desc
  compare.fun<-function(x){
    bitClade<-attach.big.matrix(desc)
    resCmp<-approximateCompare2_AVX2(bitClade@address,ntip,phy$edge,x$edge,Nnode(x),threshold=threshold+0.05)
    rm(bitClade)
    tmp<-gc()
    return(resCmp)
  }
  if(ncore>1){
    cl<-makeForkCluster(min(ncore,length(bstrees)))
    res<-parLapply(cl,bstrees,compare.fun)
    stopCluster(cl)
  }else{
    res<-lapply(bstrees,compare.fun)
  }
  unlink(clade_mat$filename)
  xorDis<-sapply(res,function(x){
    x$ref2qry[,2]
  })
  xorDis<-xorDis/clade_mat$node_size
  output<-apply(xorDis,1,function(x){
    sum(x<=threshold)/length(x)
  })
  output<-cbind(output,clade_mat$node_size)
  return(output)
}

approximate.consensus<-function(phys,prefix,threshold,ncore,CUDA_path,check=F){
  nphys<-length(phys)
  ntip<-Ntip(phys[[1]])
  res<-clade.support(phys,prefix=prefix,threshold=threshold,ncore=ncore)
  
  consensusNodes<-list() ## consensus nodes from all phys
  for(i in 1:length(res)){
    phy<-phys[[i]]
    bootstrap_res<-res[[i]]
    normailzed_dis<-bootstrap_res$dis/bootstrap_res$lib
    tmp<-apply(normailzed_dis,1,function(x){
      sum(x<=threshold)>nphys/2
    })
    if(i>1){
      for(j in 1:(i-1)){
        tmp<-tmp&(!(bootstrap_res$nodes[,j]%in%consensusNodes[[j]]&normailzed_dis[,j]<=threshold))
      }
    }
    consensusNodes[[i]]<-(1:Nnode(phy)+ntip)[tmp]
  }
  cat("No. consensus nodes: ",sapply(consensusNodes,length),fill=T)
  
  nodeRebuild<-list() ##consensus nodes from phys except for phys[[1]] need to rebuild in phys[[1]]
  for(i in 2:length(res)){
    phy<-phys[[i]]
    bootstrap_res<-res[[i]]
    normalized_dis<-bootstrap_res$dis/bootstrap_res$lib
    #retainNode<-which(bootstrap_res$nodes[,1]%in%consensusNodes[[1]]&bootstrap_res$dis[,1]<=threshold) ##nodes linked from a consensus node in phys[[1]]
    retainNode<-which(bootstrap_res$nodes[,1]%in%consensusNodes[[1]]&normalized_dis[,1]<=threshold) ##nodes linked from a consensus node in phys[[1]]
    ##if nodes link to a consensus node in phys[[1]],nodeState is the consensus node
    ##if not, nodeState is 0
    ##nodeState stores which node in phys[[1]] links the node in phy[[i]]. Therefore, there is no duplicates in nodeState
    nodeState<-c(1:ntip,rep(0,Nnode(phy)))
    nodeState[retainNode+ntip]<-bootstrap_res$nodes[retainNode,1]
    ##for each consensus node in phys[[i]],find all subclades linked from a consensus node in phys[[1]]
    ncoreRun<-min(ncore,length(consensusNodes[[i]]))
    tmpfun<-function(j){
      findMaxDes(phy$edge,j,nodeState)
    }
    if(ncoreRun>1){
      cl<-makeForkCluster(ncoreRun)
      tmp<-parLapply(cl,consensusNodes[[i]],tmpfun)
      stopCluster(cl)
    }else{
      tmp<-lapply(consensusNodes[[i]],tmpfun)
    }
    
    if(check){ ##check whether there are resulted nodes not in consensus node list of ith phy
      tmpfun<-function(x){
        sum(!x%in%c(1:ntip,retainNode+ntip))
      }
      if(ncoreRun>1){
        cl<-makeForkCluster(ncoreRun)
        aa<-parSapply(cl,tmp,tmpfun)
        stopCluster(cl)
      }else{
        aa<-sapply(tmp,tmpfun)
      }
      cat("No. descendent sets not in concensus nodes of phy",i,":",sum(aa!=0),fill=T)
    }
    
    ## change the nodes in phys[[i]] to nodes in phys[[1]]
    tmp<-lapply(tmp,function(x){
      nodeState[x]
    })
    ## joining these nodes together in phys[[1]] rebuild the consensus nodes in phys[[i]]
    nodeRebuild<-c(nodeRebuild,tmp)
    
  }
  cat("No. rebuild node sets: ",length(nodeRebuild),fill=T)
  
  if(check){ ##check whether there are rebuild nodesets not in consensus node list of 1st phy (target phy)
    cl<-makeForkCluster(ncore)
    aa<-parSapply(cl,nodeRebuild,function(x){
      sum(!x%in%c(1:ntip,consensusNodes[[1]]))
    })
    stopCluster(cl)
    cat("rebuild nodesets not in consensusNodes[[1]]:",sum(aa>0),fill=T)
  }
  
  targetNodeSize<-c(rep(1,ntip),res[[1]]$lib)
  nodeState<-c(rep(1,ntip),rep(0,Nnode(phys[[1]])))
  nodeState[consensusNodes[[1]]]<-1
  targetEdge<-collapseNode(phys[[1]]$edge,nodeState,ntip) ##mark collapsed edge by -1
  targetEdge<-targetEdge[targetEdge[,1]!=-1,] ##remove collapsed edge but do not change node idx
  
  if(check){ ##report summary of backbone phy
    cat("No. tips in phys[[1]]:", ntip,fill=T)
    cat("No. nodes in phys[[1]]:", Nnode(phys[[1]]),fill=T)
    cat("No. consensus nodes in phys[[1]]:", length(consensusNodes[[1]]),fill=T)
    cat("No. nodes+tips after collapse:",nrow(targetEdge)+1,fill=T)
  }
  
  cl<-makeForkCluster(ncore)
  nodeClean<-parLapply(cl,nodeRebuild,function(x){ 
    removeNestedNode(targetEdge,x,Nnode(phys[[1]])+ntip) ##remove nested nodes (one node is the offspring of another node)
  })
  stopCluster(cl)
  cat("No. clean node sets: ",length(nodeClean),fill=T)
  
  tmpEdge<-rebuild.nodes(
    targetEdge=targetEdge,targetNodeSize=targetNodeSize,
    ntip=ntip,nodeClean=nodeClean,
    prefix=prefix,CUDA_path=CUDA_path,
    threshold=threshold,ncore=ncore,check=F
  )
  
  ctre<-phy.wrap(tmpEdge,phys[[1]]$tip.label)
  return(ctre)
}
