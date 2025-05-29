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
#' @param weight Weights of indels.
#' @param splt The results of `split.data` from previous run.
#' @param split_idx The index of cluster to further split.
#' @param radius The minimum number of mutations between neighbors.
#' @param p_shared The threshold for the percentage of shared neighbors.
#' @param min_size The minimum size of clusters. Clusters smaller than this value will be collapsed.
#' @param batch_size The batch size during processing. Batch size has no effect on the results.
#' @param file_backed The name of backup file. If set to `NULL`, do not create backup file. If set to `""`, create an anonymous backup file.
#' @param ncore The number of threads used.
#' @param verbose Whether print detailed information.
#' 
#' 
#' @examples
#' \dontrun{
#' splt<-split.data(mat,radius=25,p_shared=0.8)
#'}
#'
#' @return Return an integer matrix containing the cluster information.

split.data<-function(mat,weight=NULL,splt=NULL,split_idx=1,
                      radius=0,p_shared=0.8,min_size=10,batch_size=1e4,file_backed=NULL,ncore=1,verbose=F){
  if(is.null(attr(mat,"bit_coded"))){
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
    # if(is.null(weight)){
    #   weight<-rep(1,attr(mat,"nloci"))
    # }
    maxChar<-attr(mat,"nloci")
  }
  rm(mat)
  
  if(is.null(splt)){
    layer<-1
    k<-1
    splt<-cbind(rep(1,nrow(binary)),rep(0,nrow(binary)))
  }else{
    layer<-ncol(splt)
    layer<-which(apply(splt,2,function(x){
      split_idx%in%x
    }))
    k<-max(splt)
    if(layer==ncol(splt)){
      splt<-cbind(splt,rep(0,nrow(binary)))
    }
  }
  
  idx<-which(splt[,layer]==split_idx)
  binary<-binary[idx,]
  if(length(idx)>batch_size){
    
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
                              p=p_shared,sqrt_nNeighbor=F,verbose=verbose)
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
    
    
  }else{
    temp<-split.mg(
      binary=binary,maxChar=maxChar,weight=weight,radius=radius,
      p_shared=p_shared,sqrt_nNeighbor=F,min_size=min_size,verbose=verbose
    )
  }
  temp[temp!=0]<-temp[temp!=0]+k
  splt[idx,layer+1]<-temp
  k<-max(splt)
  
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
    phy$edge.length[phy$edge.length==0]<-phy$edge.length[phy$edge.length==0]+1e-3
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

regularize.mat<-function(mat,ref){
  for(i in 1:ncol(mat)){
    tmp<-sort(unique(c(ref[,i],-1,0)))
    mat[,i]<-match(mat[,i],tmp)-2
  }
  return(mat)
}

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
#' @param weight Weights of indels.
#' @param noise Whether add some noise to the weights to prevent equal number of shared mutations.
#' 
#' @examples
#' \dontrun{
#' tre<-shared.mutation.join(mat)
#'}
#'
#' @return a reconstructed phylogeny.

shared.mutation.join<-function(mat,weight=NULL,noise=F){
  if(is.null(attr(mat,"bit_coded"))){
    mat<-regularize.mat(mat,mat)
    greedyfun<-greedy.partial
  }else{
    greedyfun<-greedy.partial.bit
  }
  return(greedyfun(inmat=mat,weight=weight,noise=noise))
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
#' tre<-as.ultrametric(phys)
#'}
#'
#' @return a phylogeny.

link.hcluster<-function(mg){
  output<-mg[["__M1"]]
  while(length(grep("^__M",output$tip.label))!=0){
    tip_idx<-grep("^__M",output$tip.label)[1]
    output<-bind.tree(output,mg[[output$tip.label[tip_idx]]],where=tip_idx)
  }
  return(output)
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

#' reconstruct the phylogenies using partial join algorithm
#'
#' @param mat The mutation matrix.
#' @param splt Split of data from `split.data`.
#' @param keep_cutoff The cutoff to break problematic clades.
#' @param weight Weights of indels.
#' @param noise Whether add some noise to the weights to prevent equal number of shared mutations.
#' @param bk_size If the number of OTUs is larger than this value, the function will create a backup file to save memory.
#' @param verbose Whether print detailed information.
#' 
#' @examples
#' \dontrun{
#' phys<-h.partial.join(imputed,splt,noise=F,verbose=F)
#'}
#'
#' @return a list of phylogenies.

h.partial.join<-function(mat,splt,keep_cutoff=1,weight=NULL,noise=F,bk_size=4e4,verbose=F){
  if(is.null(attr(mat,"bit_coded"))){
    mat<-regularize.mat(mat,mat)
    greedyfun<-greedy.partial
    bit_coded<-F
  }else{
    greedyfun<-greedy.partial.bit
    bit_coded<-T
  }
  splt<-cbind(splt,0)
  max_mg_idx<-max(splt)
  k<-1
  output<-list()
  for(i in ncol(splt):2-1){
    for(j in setdiff(sort(unique(splt[,i])),0)){
      idx<-which(splt[,i]==j)
      out_idx<-which(splt[,i]!=j)
      submat<-gen.mat(mat[idx,],splt[idx,i+1],bit_coded=bit_coded)
      outmat<-gen.mat(mat[out_idx,],splt[out_idx,i+1],bit_coded=bit_coded)$submat
      if(nrow(submat$submat)>bk_size){
        file_backed=""
      }else{
        file_backed=NULL
      }
      if(verbose){
        cat("Solving ",j,"th monogroup with ",nrow(submat$submat)," tips",sep="",fill=T)
      }
      partial_res<-greedyfun(
        inmat=submat$submat,
        outmat=outmat,
        keep_cutoff=keep_cutoff,
        weight=weight,
        noise=noise,
        file_backed=file_backed,
        verbose=verbose
      )
      gc()
      if(!partial_res$isbroken){
        output[[k]]<-partial_res$phys[[1]]
        names(output)[k]<-paste0("__M",j)
        k<-k+1
      }else{
        new_group<-partial_res$category[submat$idx]
        temp<-splt[idx,i+1]
        temp[new_group>0]<-new_group[new_group>0]+max_mg_idx
        splt[idx,i]<-temp
        m<-1
        while(m<=length(partial_res$phys)){
          output[[k]]<-partial_res$phys[[m]]
          names(output)[k]<-paste0("__M",m+max_mg_idx)
          k<-k+1
          m<-m+1
        }
        max_mg_idx<-max_mg_idx+max(partial_res$category)
      }
    }
  }
  return(output)
}

#' reconstruct the phylogenies using partial join algorithm in parallel.
#'
#' @param mat The mutation matrix.
#' @param splt Split of data from `split.data`.
#' @param keep_cutoff The cutoff to break problematic clades.
#' @param weight Weights of indels.
#' @param noise Whether add some noise to the weights to prevent equal number of shared mutations.
#' @param bk_size If the number of OTUs is larger than this value, the function will create a backup file to save memory.
#' @param nocre The number of threads used.
#' 
#' @examples
#' \dontrun{
#' phys<-h.partial.join.parallel(imputed,splt,noise=F,ncore=10)
#'}
#'
#' @return a list of phylogenies.

h.partial.join.parallel<-function(mat,splt,keep_cutoff=1,weight=NULL,noise=F,bk_size=4e4,ncore=2){
  if(is.null(attr(mat,"bit_coded"))){
    mat<-regularize.mat(mat,mat)
    greedyfun<-greedy.partial
    bit_coded<-F
  }else{
    greedyfun<-greedy.partial.bit
    bit_coded<-T
  }
  splt<-cbind(splt,0)
  max_mg_idx<-max(splt)
  k<-1
  output<-list()
  for(i in ncol(splt):2-1){
    j_idx<-setdiff(sort(unique(splt[,i])),0)
    cl<-makeForkCluster(ncore)
    j_idx<-sample(j_idx,length(j_idx))
    partial_res_collected<-parLapply(cl,j_idx,function(j){
      idx<-which(splt[,i]==j)
      out_idx<-which(splt[,i]!=j)
      submat<-gen.mat(mat[idx,],splt[idx,i+1],bit_coded=bit_coded)$submat
      outmat<-gen.mat(mat[out_idx,],splt[out_idx,i+1],bit_coded=bit_coded)$submat
      if(nrow(submat)>bk_size){
        file_backed=""
      }else{
        file_backed=NULL
      }
      partial_res<-greedyfun(
        inmat=submat,
        outmat=outmat,
        keep_cutoff=keep_cutoff,
        weight=weight,
        noise=noise,
        file_backed=file_backed
      )
      return(partial_res)
    })
    stopCluster(cl)
    gc()
    for(j in j_idx){
      idx<-which(splt[,i]==j)
      out_idx<-which(splt[,i]!=j)
      partial_res<-partial_res_collected[[match(j,j_idx)]]
      if(!partial_res$isbroken){
        output[[k]]<-partial_res$phys[[1]]
        names(output)[k]<-paste0("__M",j)
        k<-k+1
      }else{
        new_group<-partial_res$category[gen.mat(mat[idx,],splt[idx,i+1],return.mat=F)$idx]
        temp<-splt[idx,i+1]
        temp[new_group>0]<-new_group[new_group>0]+max_mg_idx
        splt[idx,i]<-temp
        m<-1
        while(m<=length(partial_res$phys)){
          output[[k]]<-partial_res$phys[[m]]
          names(output)[k]<-paste0("__M",m+max_mg_idx)
          k<-k+1
          m<-m+1
        }
        max_mg_idx<-max_mg_idx+max(partial_res$category)
      }
    }
  }
  return(output)
}

#' calculate the number of evolution events of the phylogeny.
#'
#' @param phy The phylogeny object.
#' @param mat The mutation matrix.
#' @param weight Weights of indels.
#' @param verbose Whether to print the number of evolution events.
#' 
#' @examples
#' \dontrun{
#' nEvo<-cal.nevo(tre,mat)
#'}
#'
#' @return an bit-coded matrix marking the evolution events in the nodes of the phylogeny.

cal.nevo<-function(phy,mat,weight=NULL,verbose=T){
  if(is.null(attr(mat,"bit_coded"))){
    mat<-rbind(mat,node.ancestral(phy,mat))
    binary<-mat2bit(mat)
    nloci<-rev(binary$loci)[1]
    temp<-diff(binary$loci)
    if(!is.null(weight)){
      weight<-unlist(lapply(which(temp!=0),function(i){
        weight[1:temp[i],i]
      }))
    }else{
      weight<-rep(1,nloci)
    }
    nEvo<-t(edgeEvo(phy$edge,binary$binary,Ntip(phy)))
  }else{
    nloci<-attr(mat,"nloci")
    mat<-rbind(mat,node.ancestral.bit(phy,mat))
    if(is.null(weight)){
      weight<-rep(1,nloci)
    }
    nEvo<-apply(phy$edge,1,function(x){
      locusXOR(mat[x[1],],mat[x[2],])
    })
    nEvo<-nEvo[,match(1:(Ntip(phy)+Nnode(phy)),phy$edge[,2])]
    nEvo[,Ntip(phy)+1]<-mat[Ntip(phy)+1,]
  }
  if(verbose){
    total_nEvo<-sum(apply(nEvo,2,locusSUM,weight=weight))
    cat("No. evolution:",total_nEvo,fill=T)
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

#' add the length of edges to the phylogeny according the number of evolution events and collapse nodes without evolution.
#'
#' @param phy The phylogeny object.
#' @param nEvo The results from `cal.nevo`.
#' @param collapse Whether to collapse the nodes. If set to False, the function just add the length of edges to the phylogeny.
#' 
#' @examples
#' \dontrun{
#' tre<-collapse.infoless(tre,nEvo)
#'}
#'
#' @return a phylogeny.

collapse.infoless<-function(phy,nEvo,collapse=T){
  phy<-add.nEvo(phy,nEvo)
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

find.dup<-function(mat,ncore=1){
  if(is.null(attr(mat,"bit_coded"))){
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

# iter.impute.packed2<-function(mat,weight=NULL,n_neighbor=10:2,p_mutation=0.5,n_iter=100,cutoff=max(floor(nrow(mat)/100),100),verbose=F,ncore=1){
#   binary<-mat2bit(mat)
#   nloci<-rev(binary$loci)[1]
#   
#   if (!is.null(weight)&class(weight)[1]=="matrix") {
#     temp <- diff(binary$loci)
#     weight <- unlist(lapply(which(temp != 0), function(i) {
#       weight[1:temp[i], i]
#     }))
#   }
#   
#   id<-rownames(mat)
#   maxChar<-ncol(mat)
#   
#   n_mask<-sum(mat==-1)
#   total_imputed<-0
#   k<-0
#   neighbor_idx<-1
#   initial<-T
#   reinitial<-0
#   
#   if(is.null(weight)){
#     fun1<-function(i){
#       findNeighborUnweightedParallel(i=i-1,mat=binary,maxChar=maxChar,
#                                      n=n_neighbor[neighbor_idx],p=p_mutation)
#     }
#     fun2<-function(i){
#       findNeighborRangeUnweighted(i=i-1,mat=binary,maxChar=maxChar,
#                                   range=range_found,offset=0,
#                                   n=n_neighbor[neighbor_idx],p=p_mutation)
#     }
#   }else{
#     fun1<-function(i){
#       findNeighborWeightedParallel(i=i-1,mat=binary,maxChar=maxChar,
#                                    weight=weight,n=n_neighbor[neighbor_idx],p=p_mutation)
#     }
#     fun2<-function(i){
#       findNeighborRangeWeighted(i=i-1,mat=binary,maxChar=maxChar,
#                                 range=range_found,offset=0,
#                                 weight=weight,n=n_neighbor[neighbor_idx],p=p_mutation)
#     }
#   }
#   impute.fun1<-function(i){
#     imputeMat(i-1,mat,range_found,0)
#   }
#   impute.fun2<-function(i){
#     imputeMat(i-1,mat,neighbor_found,0)
#   }  
#   {
#     if(verbose){
#       cat(as.character(Sys.time()),": iteration ",k+1,", No. neighbor ",n_neighbor[neighbor_idx],fill=T)
#     }
#     binary<-mat2bit(mat)$binary
#     if(ncore>1){
#       cl<-makeForkCluster(ncore)
#       neighbor_found<-t(parSapply(cl,1:nrow(binary),fun1))
#       mat<-t(parSapply(cl,1:nrow(mat),impute.fun2))
#       stopCluster(cl)
#     }else{
#       neighbor_found<-t(sapply(1:nrow(binary),fun1))
#       mat<-t(sapply(1:nrow(mat),impute.fun2))
#     }
#     k<-k+1
#     temp<-n_mask-sum(mat==-1)
#     total_imputed<-total_imputed+temp
#     if(verbose){
#       cat("Impute ",temp," entries.",fill=T)
#     }
#     n_mask<-n_mask-temp
#   }
#   
#   while(k<n_iter&neighbor_idx<=length(n_neighbor)){
#     if(verbose){
#       cat(as.character(Sys.time()),": iteration ",k+1,", No. neighbor ",n_neighbor[neighbor_idx],fill=T)
#     }
#     binary<-mat2bit(mat)$binary
#     if(initial){
#       if(verbose){
#         cat("initializing",fill=T)
#       }
#       if(ncore>1){
#         cl<-makeForkCluster(ncore)
#         range_found<-t(parSapply(cl,1:nrow(binary),fun1))
#         mat<-t(parSapply(cl,1:nrow(mat),impute.fun1))
#         stopCluster(cl)
#       }else{
#         range_found<-t(sapply(1:nrow(binary),fun1))
#         mat<-t(sapply(1:nrow(mat),impute.fun1))
#       }
#       initial<-F
#     }else{
#       if(ncore>1){
#         cl<-makeForkCluster(ncore)
#         neighbor_found<-t(parSapply(cl,1:nrow(mat),fun2))
#         mat<-t(parSapply(cl,1:nrow(mat),impute.fun2))
#         stopCluster(cl)
#       }else{
#         neighbor_found<-t(sapply(1:nrow(mat),fun2))
#         mat<-t(sapply(1:nrow(mat),impute.fun2))
#       }
#     }
#     k<-k+1
#     temp<-n_mask-sum(mat==-1)
#     total_imputed<-total_imputed+temp
#     if(verbose){
#       cat("Impute ",temp," entries.",fill=T)
#     }
#     n_mask<-n_mask-temp
#     if(temp<cutoff){
#       neighbor_idx<-neighbor_idx+1
#       reinitial<-reinitial+1
#       if(reinitial>5){
#         initial<-T
#         reinitial<-0
#       }
#     }
#   }
#   rownames(mat)<-id
#   if(verbose){
#     cat("Impute ",k," times with ",total_imputed," entries in total.",sep="",fill=T)
#   }
#   return(mat)
# }

write.mat2bin<-function(mat,file){
  output<-rev(dim(mat))
  output<-c(output,mat)
  writeBin(object=output,con=file)
}

read.bin2mat<-function(file){
  d<-readBin(con=file,what="integer",n=2)
  mat<-readBin(con=file,what="integer",n=d[1]*d[2]+2)
  mat<-matrix(mat[c(-1,-2)],d[2],d[1])
  return(mat)
}
