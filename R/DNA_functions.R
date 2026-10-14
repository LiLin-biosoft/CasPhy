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
    # tmp<-table(rootVote[retain,1])/length(retain)
    tmp<-tapply(rootVote[retain,3],rootVote[retain,1],sum)/sum(rootVote[retain,3])
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
  
  ## relabel node idx of phy2
  transition<-rep(0,Nnode(phy2))
  transitionIdx<-unique(phy2$edge[phy2$edge>Ntip(phy2)])-Ntip(phy2)
  transition[transitionIdx]<-1:Nnode(phy2)+Ntip(phy2) ##old idx to new idx
  phy2$node.label[transition[1:Nnode(phy2)]-Ntip(phy2)]<-phy2$node.label
  phy2$edge[phy2$edge>Ntip(phy2)]<-transition[phy2$edge[phy2$edge>Ntip(phy2)]-Ntip(phy2)]
  if(length(break_node)==0){
    return(list(phy=phy2,break_edge=matrix(,0,2),break_node=integer()))
  }
  ## relabel node idx of phy2
  
  ##tranlate node idx in phy to node idx in phy2
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
  ##tranlate node idx in phy to node idx in phy2
  return(list(phy=phy2,break_edge=break_edge2,break_node=break_node2))
}

partial.join.vft<-function(j,mat_desc,tipname,matfile,splt,layer,prefix,cuda_path,sample_size=0.1,
                           pct_threshold=0.001,nbp_threshold=200,
                           vft_cmd="VeryFastTree -nosupport -fastest -threads 8",
                           outgroup_mode=0,break_mode=1){
  ## outgroup_mode=0 use true OTUs
  ## outgroup_mode=1 use pseudo OTUs in the same layer
  ## outgroup_mode=2 use pseudo OTUs weighted by the size
  
  ## break_mode=0 break primary root
  ## break_mode=1 break all candidate root
  ## break_mode=2 break primary root and random edge
  
  mat<-attach.big.matrix(mat_desc)
  inIdx<-which(splt[,layer]==j)
  outIdx<-which(splt[,layer]!=j)
  if(outgroup_mode==1|outgroup_mode==2){
    nodename<-paste0("__M",setdiff(unique(splt[outIdx,layer+1]),0))
    outIdx<-c(
      outIdx[splt[outIdx,layer+1]==0],
      match(nodename,tipname)
    )
  }
  
  inmat<-gen.mat.big(mat,tipname,splt[,layer+1],inIdx=inIdx)
  tip2pseudo<-inmat$idx
  inmat<-inmat$submat
  tipLabel<-colnames(inmat)
  
  ##run VeryFastTree
  write.mat2fa(inmat,file=paste0(prefix,"__M",j,".fasta"),randomize_ambiguous=F)
  vft_tmp_file<-c(
    paste0(prefix,"__M",j,".fasta",sep=""),
    paste0(prefix,"__M",j,".nwk",sep=""),
    paste0(prefix,"__M",j,".vftlog")
  )
  system(paste0(
    vft_cmd," ",
    "-nt ",prefix,"__M",j,".fasta ",
    " 1>",prefix,"__M",j,".nwk 2>",prefix,"__M",j,".vftlog"
  ))
  phy<-read.tree(paste0(prefix,"__M",j,".nwk"))
  phy<-order.tip(phy,tipLabel)
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-cbind(inmat,node_mat)
  
  if(length(outIdx)==0){
    unlink(vft_tmp_file)
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
    
    voteRoot_tmp_file<-c(
      paste0(prefix,"__M",j,".nodeState"),
      paste0(prefix,"__M",j,".outmatIdx"),
      paste0(prefix,"__M",j,".edge"),
      paste0(prefix,"__M",j,".rootVote")
    )
    
    system(paste0(cuda_path,"voteRoot2 ",matfile," ",prefix,"__M",j))
    rootVote<-read.bin2mat(paste0(prefix,"__M",j,".rootVote"))
    rootVote[,1]<-rootVote[,1]+1L
    unlink(voteRoot_tmp_file)
  }
  if(outgroup_mode==0|outgroup_mode==1){
    rootVote<-cbind(rootVote,1)
  }else{
    tmp<-table(splt[,layer+1])
    names(tmp)<-paste0("__M",names(tmp))
    ntip<-nrow(splt)
    outSize<-c(
      rep(1,sum(outIdx<=ntip)),
      tmp[tipname[outIdx[outIdx>ntip]]]
    )
    rootVote<-cbind(rootVote,outSize)
  }
  
  breakRes<-get.break.node(rootVote,cutoff=nbp_threshold,pct_threshold=pct_threshold)
  root_node<-breakRes$root_node
  break_node<-breakRes$break_node
  rootedRes<-break.pp(phy,root_node,break_node)
  phy<-rootedRes$phy
  break_node<-rootedRes$break_node
  if(break_mode==0){
    break_node<-integer()
  }else if(break_mode==2){
    break_node<-sample(2:Nnode(phy)+Ntip(phy),length(break_node))
  }
  nodes_retain<-breakUnrootedPhy(phy$edge,Ntip(phy),Nnode(phy),breakNode=break_node)
  nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
  if(length(nodes_retain)==0){
    unlink(vft_tmp_file)
    return(list(phys=list(),category=rep(0,length(inIdx)),ansState=NULL,isbroken=T,inIdx=inIdx))
  }
  t2<-Sys.time()
  
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-node_mat[,nodes_retain-Ntip(phy)]
  tmp<-phy.break(phy,nodes_retain)
  broken_phys<-tmp$phys
  category<-tmp$category
  unlink(vft_tmp_file)
  return(list(phys=broken_phys,category=category[tip2pseudo],ansState=node_mat,isbroken=T,inIdx=inIdx))
}

h.partial.join.vft<-function(input_file,splt,prefix,sample_size=1,pct_threshold=0.001,nbp_threshold=NULL,file_backed=F,
                             vft_cmd="VeryFastTree -nosupport -fastest -threads 8",cuda_path,ncore){
  
  cat(as.character(Sys.time()),"pvft start",fill=T)
  ntip<-nrow(splt)
  mat<-read.bin2mat(input_file)
  if(is.null(nbp_threshold)){
    nbp_threshold<-ceiling(nrow(mat)*4*0.05)
  }
  
  matfile<-paste0(prefix,".bin.tmp")
  write.mat2bin(mat,file=matfile,extraCols=floor(ntip*0.5))
  rm(mat)
  gc()
  
  if(!all(splt[,ncol(splt)]==0)){
    splt<-cbind(splt,0)
  }
  
  if(file_backed){
    fpath<-gsub("[^/]+$","",prefix)
    fname<-paste0(gsub(fpath,"",prefix),".bk")
    mat<-read.bin2bigmat(
      matfile,
      backingfile=fname,
      backingpath=fpath,
      descriptorfile=paste0(fname,".desc")
    )
  }else{
    mat<-read.bin2bigmat(matfile)
  }
  mat_desc<-describe(mat)
  tipname<-paste0("C",1:ntip)
  phys<-list()
  for(layer in ncol(splt):3-1){
    clustIdx<-setdiff(sort(unique(splt[,layer])),0)
    partial_res_collected<-list()
    nProcessed<-0
    clusterSize<-sapply(clustIdx,function(x){
      inTip<-splt[splt[,layer]==x,layer+1]
      return(get.phySize(inTip))
    })
    clustIdx<-clustIdx[order(clusterSize,decreasing=T)]
    
    cat(
      as.character(Sys.time()),
      "layer:",layer,
      "No. clusters:",length(clustIdx),
      "max size:",max(clusterSize),"\n",fill=F
    )
    ncoreIdx<-1
    while(nProcessed<length(clustIdx)){
      if(ncoreIdx<length(ncore)){
        ncoreRun<-ncore[ncoreIdx]
        ncoreIdx<-ncoreIdx+1
      }else{
        ncoreRun<-ncore[ncoreIdx]
      }
      
      cl<-makeForkCluster(ncoreRun)
      tmp_res<-parLapply(cl,clustIdx[(nProcessed+1):min((nProcessed+ncoreRun),length(clustIdx))],partial.join.vft,
                         mat_desc=mat_desc,
                         tipname=tipname,
                         matfile=matfile,
                         splt=splt,
                         layer=layer,
                         prefix=prefix,
                         cuda_path=cuda_path,
                         sample_size=sample_size,
                         pct_threshold=pct_threshold,
                         nbp_threshold=nbp_threshold,
                         vft_cmd=vft_cmd,
                         outgroup_mode=0
      )
      stopCluster(cl)
      partial_res_collected<-c(partial_res_collected,tmp_res)
      nProcessed<-min((nProcessed+ncoreRun),length(clustIdx))
      rm(tmp_res)
      tmp<-gc()
      cat(as.character(Sys.time()),"No. clusters remain:",length(clustIdx)-nProcessed,fill=T)
    }
    
    cat(as.character(Sys.time()),"summerizing layer ",layer,fill=T)
    res_summarized<-layer.summarize(partial_res_collected,clustIdx,splt,layer)
    if(ncol(res_summarized$ansMat)+length(tipname)>ncol(mat)){
      cat(as.character(Sys.time()),"ans cache filled",fill=T)
      return(0)
    }
    phys<-c(phys,res_summarized$phys)
    mat[,1:ncol(res_summarized$ansMat)+length(tipname)]<-res_summarized$ansMat
    write.mat2bin(as.vector(res_summarized$ansMat),matfile,append=length(tipname)*nrow(mat))
    tipname<-c(tipname,colnames(res_summarized$ansMat))
    cat("ans cache used:",length(tipname)-ntip,fill=T)
    splt<-res_summarized$splt  
    cat(as.character(Sys.time())," layer ",layer," finished",fill=T)
  }
  
  ##########################################
  phySize<-get.phySize(splt[,2])
  cat(as.character(Sys.time()),"solving backbone. No. tips:",phySize,fill=T)
  partial_res_collected<-list(
    partial.join.vft(
      j=1,
      mat_desc=mat_desc,
      tipname=tipname,
      matfile=matfile,
      splt=splt,
      layer=1,
      prefix=prefix,
      cuda_path=cuda_path,
      sample_size=sample_size,
      pct_threshold=pct_threshold,
      nbp_threshold=nbp_threshold,
      vft_cmd=gsub("-threads\\s+[0-9]+",paste0("-threads ",rev(ncore)[1]),vft_cmd)
    )
  )
  res_summarized<-layer.summarize(
    partial_res_collected=partial_res_collected,
    clustIdx=1,
    splt=splt,
    layer=1
  )
  phys<-c(phys,res_summarized$phys)
  cat(as.character(Sys.time()),"build tree finished.",fill=T)
  rm(mat)
  gc()
  if(file_backed){
    unlink(paste0(prefix,".bk"))
    unlink(paste0(prefix,".bk.desc"))
  }
  unlink(matfile)
  phys_bak<-phys
  return(phys)
  
  cat(as.character(Sys.time()),"solving layer 1",fill=T)
  nRound<-0
  phySize<-get.phySize(splt[,2])
  # while(phySize>sizeBound){
  while(nRound<=0){
    nRound<-nRound+1
    cat(as.character(Sys.time()),"No. round:",nRound,"No. tips:",phySize,fill=T)
    matTmp<-gen.mat.big(mat,tipname,splt[,2],1:nrow(splt))
    idxTmp<-matTmp$idx
    res<-autoSplit(
      mat=matTmp$submat,prefix=prefix,radius_initial=radius_initial-layerStep*(nRound-1),
      fixed=F,step=tryStep,radius_bound=tryStep,
      sizeBound=c(100000L,0L),cuda_path=cuda_path
    )
    rm(matTmp)
    tmp<-gc()
    spltTmp<-res[[2]]
    spltTmp<-regularize.split(spltTmp,min_size=10)
    nPseudoTip<-get.phySize(spltTmp)
    if(nPseudoTip<phySize*0.5){
      spltTmp<-addOffset(spltTmp,max(splt))
      spltTmp<-spltTmp[idxTmp]
      splt<-cbind(1,spltTmp,splt[,2:ncol(splt)])
    }else{
      cat(as.character(Sys.time())," spliting failed in round:",nRound,"No. tips:",phySize,fill=T)
      break
    }
    ############################
    clustIdx<-setdiff(sort(unique(splt[,2])),0)
    clusterSize<-sapply(clustIdx,function(x){
      inTip<-splt[splt[,2]==x,3]
      return(get.phySize(inTip))
    })
    clustIdx<-clustIdx[order(clusterSize,decreasing=T)]
    cat(as.character(Sys.time()),
        "round:",nRound,
        "radius:",res[[1]],
        "No. clusters:",length(clustIdx),
        "max size:",max(clusterSize),"\n",fill=F
    )
    partial_res_collected<-list()
    nProcessed<-0
    while(nProcessed<length(clustIdx)){
      if(nProcessed==0){
        ncoreRun<-floor(ncore/2)
      }else{
        ncoreRun<-ncore
      }
      cl<-makeForkCluster(ncoreRun)
      tmp_res<-parLapply(cl,clustIdx[(nProcessed+1):min((nProcessed+ncoreRun),length(clustIdx))],partial.join.vft,
                         mat_desc=mat_desc,
                         tipname=tipname,
                         matfile=matfile,
                         splt=splt,
                         layer=2,
                         prefix=prefix,
                         cuda_path=cuda_path,
                         sample_size=sample_size,
                         pct_threshold=pct_threshold,
                         nbp_threshold=nbp_threshold,
                         outgroup_mode=0
      )
      stopCluster(cl)
      partial_res_collected<-c(partial_res_collected,tmp_res)
      nProcessed<-min((nProcessed+ncoreRun),length(clustIdx))
      rm(tmp_res)
      tmp<-gc()
      cat(as.character(Sys.time()),"No. clusters remain:",length(clustIdx)-nProcessed,fill=T)
    }
    cat(as.character(Sys.time()),"summerizing round ",nRound,fill=T)
    res_summarized<-layer.summarize(partial_res_collected,clustIdx,splt,2)
    if(ncol(res_summarized$ansMat)+length(tipname)>ncol(mat)){
      cat(as.character(Sys.time()),"ans cache filled",fill=T)
      return(0)
    }
    phys<-c(phys,res_summarized$phys)
    mat[,1:ncol(res_summarized$ansMat)+length(tipname)]<-res_summarized$ansMat
    write.mat2bin(as.vector(res_summarized$ansMat),matfile,append=length(tipname)*nrow(mat))
    tipname<-c(tipname,colnames(res_summarized$ansMat))
    cat("ans cache used:",length(tipname)-ntip,fill=T)
    splt<-res_summarized$splt
    phySize<-get.phySize(splt[,2])
  }
  cat(as.character(Sys.time()),"building final tree with",phySize,"tips.",fill=T)
  partial_res_collected<-list(
    partial.join.vft(
      j=1,
      mat_desc=mat_desc,
      tipname=tipname,
      matfile=matfile,
      splt=splt,
      layer=1,
      prefix=prefix,
      cuda_path=cuda_path,
      sample_size=sample_size,
      pct_threshold=pct_threshold,
      nbp_threshold=nbp_threshold,
      vft_cmd=gsub("-threads\\s+[0-9]+",paste0("-threads ",rev(ncore)[1]),vft_cmd)
    )
  )
  res_summarized<-layer.summarize(
    partial_res_collected=partial_res_collected,
    clustIdx=1,
    splt=splt,
    layer=1
  )
  phys<-c(phys,res_summarized$phys)
  return(list(phys_bak,phys))
}

epa2break<-function(epa_res,outSize){
  like_weight_ratio_sum<-tapply(
    epa_res$placement_positions$like_weight_ratio*outSize[epa_res$placement_positions$placement_id],
    epa_res$placement_positions$edge_num,
    sum
  )
  #  insert_edge<-epa_res$edge_key[1,match(as.integer(names(like_weight_ratio_sum)),epa_res$edge_key[2,])]
  insert_edge<-as.integer(names(like_weight_ratio_sum))
  insert_edge<-insert_edge[order(like_weight_ratio_sum,decreasing=T)]
  insert_node<-epa_res$arbre$edge[insert_edge,2]
  insert_node<-insert_node[1:ceiling(length(insert_node*0.1))]
  root_node<-insert_node[1]
  break_node<-setdiff(insert_node,root_node)
  return(list(root_node=root_node,break_node=break_node))
}

partial.join.raxml<-function(j,mat_desc,tipname,matfile,splt,layer,prefix,sample_size=0.1,
                             raxml_cmd="raxml-ng --redo --model GTR+G --force perf_threads --tree pars{1} --threads auto {4}",
                             epa_cmd="epa-ng --redo --threads 4",
                             outgroup_mode=0,break_mode=1){
  mat<-attach.big.matrix(mat_desc)
  inIdx<-which(splt[,layer]==j)
  outIdx<-which(splt[,layer]!=j)
  if(outgroup_mode==1|outgroup_mode==2){
    nodename<-paste0("__M",setdiff(unique(splt[outIdx,layer+1]),0))
    outIdx<-c(
      outIdx[splt[outIdx,layer+1]==0],
      match(nodename,tipname)
    )
  }
  
  inmat<-gen.mat.big(mat,tipname,splt[,layer+1],inIdx=inIdx)
  tip2pseudo<-inmat$idx
  inmat<-inmat$submat
  tipLabel<-colnames(inmat)
  
  ##run RAxML
  write.mat2fa(inmat,file=paste0(prefix,"__M",j,".fasta"),randomize_ambiguous=F)
  raxml_tmp_file<-c(
    paste0(prefix,"__M",j,".fasta",sep=""),
    paste0(prefix,"__M",j,".raxml.rba"),
    paste0(prefix,"__M",j,".raxml.reduced.phy"),
    paste0(prefix,"__M",j,".raxml.bestTreeCollapsed"),
    paste0(prefix,"__M",j,".raxml.bestTree"),
    paste0(prefix,"__M",j,".raxml.startTree"),
    paste0(prefix,"__M",j,".raxml.bestModel"),
    paste0(prefix,"__M",j,".raxml.log"),
    paste0(prefix,"__M",j,".raxlog")
  )
  system(paste0(
    raxml_cmd," ",
    "--msa ",prefix,"__M",j,".fasta ",
    "--prefix ",prefix,"__M",j,
    " >",prefix,"__M",j,".raxlog"
  ))
  phy<-read.tree(paste0(prefix,"__M",j,".raxml.bestTree"))
  phy<-order.tip(phy,tipLabel)
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-cbind(inmat,node_mat)
  t2<-Sys.time()
  
  if(length(outIdx)==0){
    unlink(raxml_tmp_file)
    return(list(phys=list(phy),category=rep(1,length(inIdx)),ansState=node_mat[,1],isbroken=F))
  }
  
  nSample<-max(16384,floor(length(outIdx)*sample_size))
  if(length(outIdx)>nSample){
    outIdx<-sort(sample(outIdx,nSample))
  }
  
  ##run EPA
  if(outgroup_mode==0|outgroup_mode==1){
    outSize<-rep(1,length(outIdx))
  }else{
    tmp<-table(splt[,layer+1])
    names(tmp)<-paste0("__M",names(tmp))
    ntip<-nrow(splt)
    outSize<-c(
      rep(1,sum(outIdx<=ntip)),
      tmp[tipname[outIdx[outIdx>ntip]]]
    )
  }
  
  outmat<-mat[,outIdx]
  colnames(outmat)<-tipname[outIdx]
  write.mat2fa(outmat,file=paste0(prefix,"__M",j,"_out.fasta"),randomize_ambiguous=F)
  epa_tmp_file<-c(
    paste0(prefix,"__M",j,"_epa"),
    paste0(prefix,"__M",j,"_out.fasta"),
    paste0(prefix,"__M",j,".epalog")
  )
  unlink(epa_tmp_file[1],recursive=T)
  dir.create(epa_tmp_file[1])
  system(paste0(
    epa_cmd," ",
    "--ref-msa ",prefix,"__M",j,".fasta ",
    "--tree ",prefix,"__M",j,".raxml.bestTree ",
    "--query ",prefix,"__M",j,"_out.fasta ",
    "--model ",prefix,"__M",j,".raxml.bestModel ",
    "--out-dir ",prefix,"__M",j,"_epa",
    " >",prefix,"__M",j,".epalog"
  ))
  epa_res<-read_jplace(paste0(prefix,"__M",j,"_epa/epa_result.jplace"))
  #unlink(paste0(prefix,"__M",j,"_epa"),recursive=T)
  epa_res$arbre<-order.tip(epa_res$arbre,tipLabel)
  breakRes<-epa2break(epa_res,outSize)
  rootedRes<-break.pp(
    epa_res$arbre,
    breakRes$root_node,
    breakRes$break_node
  )
  unlink(raxml_tmp_file)
  unlink(epa_tmp_file,recursive=T)
  phy<-rootedRes$phy
  break_node<-rootedRes$break_node
  if(break_mode==0){
    break_node<-integer()
  }else if(break_mode==2){
    break_node<-sample(2:Nnode(phy)+Ntip(phy),length(break_node))
  }
  nodes_retain<-breakUnrootedPhy(phy$edge,Ntip(phy),Nnode(phy),breakNode=break_node)
  nodes_retain<-sort(nodes_retain[nodes_retain>Ntip(phy)])
  if(length(nodes_retain)==0){
    return(list(phys=list(),category=rep(0,length(inIdx)),ansState=NULL,isbroken=T,inIdx=inIdx))
  }
  
  node_mat<-nucleoAncestral(inmat,phy$edge,Ntip(phy),Nnode(phy))$nodeState
  node_mat<-node_mat[,nodes_retain-Ntip(phy)]
  tmp<-phy.break(phy,nodes_retain)
  broken_phys<-tmp$phys
  category<-tmp$category
  return(list(phys=broken_phys,category=category[tip2pseudo],ansState=node_mat,isbroken=T,inIdx=inIdx))
}

h.partial.join.raxml<-function(input_file,splt,prefix,sample_size=1,
                               raxml_cmd="raxml-ng --redo --model GTR+G --force perf_threads --tree pars{1} --threads auto {4}",
                               epa_cmd="epa-ng --redo --threads 4",
                               ncore){
  cat(as.character(Sys.time()),"prax start",fill=T)
  ntip<-nrow(splt)
  mat<-read.bin2mat(input_file)
  matfile<-paste0(prefix,".bin.tmp")
  write.mat2bin(mat,file=matfile,extraCols=floor(ntip*0.5))
  rm(mat)
  gc()
  
  if(!all(splt[,ncol(splt)]==0)){
    splt<-cbind(splt,0)
  }
  
  mat<-read.bin2bigmat(matfile)
  mat_desc<-describe(mat)
  tipname<-paste0("C",1:ntip)
  phys<-list()
  for(layer in ncol(splt):3-1){
    clustIdx<-setdiff(sort(unique(splt[,layer])),0)
    partial_res_collected<-list()
    nProcessed<-0
    clusterSize<-sapply(clustIdx,function(x){
      inTip<-splt[splt[,layer]==x,layer+1]
      return(get.phySize(inTip))
    })
    clustIdx<-clustIdx[order(clusterSize,decreasing=T)]
    
    cat(
      as.character(Sys.time()),
      "layer:",layer,
      "No. clusters:",length(clustIdx),
      "max size:",max(clusterSize),"\n",fill=F
    )
    while(nProcessed<length(clustIdx)){
      ncoreRun<-ncore
      if(ncoreRun>1){
        cl<-makeForkCluster(ncoreRun)
        tmp_res<-parLapply(cl,clustIdx[(nProcessed+1):min((nProcessed+ncoreRun),length(clustIdx))],partial.join.raxml,
                           mat_desc=mat_desc,
                           tipname=tipname,
                           matfile=matfile,
                           splt=splt,
                           layer=layer,
                           prefix=prefix,
                           sample_size=sample_size,
                           raxml_cmd=raxml_cmd,
                           epa_cmd=epa_cmd
        )
        stopCluster(cl)
      }else{
        tmp_res<-list(partial.join.raxml(
          j=clustIdx[nProcessed+1],
          mat_desc=mat_desc,
          tipname=tipname,
          matfile=matfile,
          splt=splt,
          layer=layer,
          prefix=prefix,
          sample_size=sample_size,
          raxml_cmd=raxml_cmd,
          epa_cmd=epa_cmd
        ))
      }
      partial_res_collected<-c(partial_res_collected,tmp_res)
      nProcessed<-min((nProcessed+ncoreRun),length(clustIdx))
      rm(tmp_res)
      tmp<-gc()
      cat(as.character(Sys.time()),"No. clusters remain:",length(clustIdx)-nProcessed,fill=T)
    }
    
    cat(as.character(Sys.time()),"summerizing layer ",layer,fill=T)
    res_summarized<-layer.summarize(partial_res_collected,clustIdx,splt,layer)
    if(ncol(res_summarized$ansMat)+length(tipname)>ncol(mat)){
      cat(as.character(Sys.time()),"ans cache filled",fill=T)
      return(0)
    }
    phys<-c(phys,res_summarized$phys)
    mat[,1:ncol(res_summarized$ansMat)+length(tipname)]<-res_summarized$ansMat
    write.mat2bin(as.vector(res_summarized$ansMat),matfile,append=length(tipname)*nrow(mat))
    tipname<-c(tipname,colnames(res_summarized$ansMat))
    cat("ans cache used:",length(tipname)-ntip,fill=T)
    splt<-res_summarized$splt  
    cat(as.character(Sys.time())," layer ",layer," finished",fill=T)
  }
  
  ##########################################
  
  phySize<-get.phySize(splt[,2])
  cat(as.character(Sys.time()),"solving backbone. No. tips:",phySize,fill=T)
  partial_res_collected<-list(
    partial.join.raxml(
      j=1,
      mat_desc=mat_desc,
      tipname=tipname,
      matfile=matfile,
      splt=splt,
      layer=1,
      prefix=prefix,
      sample_size=sample_size,
      raxml_cmd=gsub("--threads auto \\{[0-9]+\\}",paste0("--threads auto {",ncore,"}"),raxml_cmd),
      epa_cmd=epa_cmd
    )
  )
  res_summarized<-layer.summarize(
    partial_res_collected=partial_res_collected,
    clustIdx=1,
    splt=splt,
    layer=1
  )
  phys<-c(phys,res_summarized$phys)
  unlink(matfile)
  return(phys)
}

read_jplace<-function (jplace_file, full = TRUE){
  json_data <- fromJSON(jplace_file)
  arbre <- json_data$tree
  arbre2 <- gsub("}", "", gsub("{", "#", 
                               arbre, fixed = TRUE), fixed = TRUE)
  arbre3 <- read.tree(text = gsub("#", "", gsub(":[0-9].[0-9]+#", 
                                                ":", gsub("e-", "", arbre2))))
  edge_key <- rbind(1:nrow(arbre3$edge), arbre3$edge.l)
  arbre <- gsub("}", "#", gsub("{", "#", 
                               arbre, fixed = TRUE), fixed = TRUE)
  arbre <- read.tree(text = gsub("#[0-9]+#", "", 
                                 arbre))
  out <- list(arbre, edge_key, json_data$tree)
  names(out) <- c("arbre", "edge_key", "original_tree")
  if (full) {
    a <- json_data$placements[, 1]
    for (i in 1:length(a)) {
      a[[i]] <- cbind(rep(i, nrow(a[[i]])), a[[i]])
    }
    a <- as.data.frame(do.call(rbind, a), stringsAsFactors = FALSE)
    colnames(a) <- c("placement_id", json_data$fields)
    col_num = c(1, 3:ncol(a))
    a[, col_num] <- apply(a[, col_num], 2, as.numeric)
    out$placement_positions <- a
    b <- json_data$placements[, 2]
    if (class(b[[1]]) == "matrix") {
      for (i in 1:length(b)) {
        b[[i]] <- cbind(rep(i, nrow(b[[i]])), b[[i]])
      }
      b2 <- as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
      colnames(b2) <- c("placement_id", "name", 
                        "nm")
      b2[, 1] <- as.numeric(b2[, 1])
      b2[, 3] <- as.numeric(b2[, 3])
    }
    if (class(b[[1]]) != "matrix") {
      for (i in 1:length(b)) {
        b[[i]] <- cbind(rep(i, length(b[[i]])), b[[i]])
      }
      b2 <- as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
      colnames(b2) <- c("placement_id", "name")
      b2[, 1] <- as.numeric(b2[, 1])
    }
    out$multiclass <- b2
    pplacer_branch_id <- out$placement_positions$edge_num
    out$placement_positions$edge_num <- out$edge_key[1, match(pplacer_branch_id, 
                                                              out$edge_key[2, ])]
    if (ncol(out$placement_positions) == 7) 
      out$placement_positions <- out$placement_positions[, 
                                                         c(1, 4, 5, 6, 3, 7, 2)]
    if (ncol(out$placement_positions) == 6) 
      out$placement_positions <- out$placement_positions[, 
                                                         c(1, 2, 4, 3, 5, 6)]
    out$run <- as.character(json_data$metadata$invocation[1])
    class(out) <- "jplace"
  }
  out
}