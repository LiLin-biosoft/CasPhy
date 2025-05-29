#include <Rcpp.h>

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

using namespace Rcpp;

// [[Rcpp::export]]
int mutationNested(IntegerVector des,IntegerVector ans){
  int i;
  int output=0;
  for(i=0;i<des.size();i++){
    if(ans(i)>0&&des(i)>=0&&ans(i)!=des(i)){
      output++;
    }
  }
  return output;
}

// [[Rcpp::export]]
int mutationNestedBit(IntegerVector des,IntegerVector ans){
  int i,j;
  int output=0;
  for(i=0;i<des.size();i++){
    int tmp=ans(i)&(~des(i));
    for(j=0;j<32;j++){
      if(((tmp>>j)&1)==1){
        output++;
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
int transitionRoot(IntegerVector transition,int idx){
  if(idx==-1){
    return -1;
  }else{
    int root=idx;
    while(transition(root)!=root&&transition(root)!=-1){
      root=transition(root);
    }
    return root;
  }
}

// [[Rcpp::export]]
IntegerVector matCoding(IntegerMatrix mat,IntegerVector loci){
  int i,j;
  int nbin=((loci(loci.size()-1)-1)>>5)+1;
  IntegerMatrix output(mat.nrow(),nbin);
  for(i=0;i<mat.nrow();i++){
    for(j=0;j<mat.ncol();j++){
      if(mat(i,j)>0){
        int bin_x=(loci(j)+mat(i,j)-1)>>5;
        int bin_y=(loci(j)+mat(i,j)-1)&31;
        output(i,bin_x)=output(i,bin_x)|(1<<bin_y);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
LogicalVector int2binary(IntegerVector ve,int output_len=0){
  int i,j;
  int len=ve.size()*32;
  if(output_len!=0){
    len=output_len;
  }
  LogicalVector output(len);
  for(i=0;i<ve.size();i++){
    for(j=0;j<32&&(i*32+j)<len;j++){
      if(((ve(i)>>j)&1)==1){
        output(i*32+j)=true;
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> int2idx(IntegerVector ve){
  int i,j;
  std::vector<int> output;
  for(i=0;i<ve.size();i++){
    for(j=0;j<32;j++){
      if(((ve(i)>>j)&1)==1){
        output.push_back(i*32+j);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector idx2int(IntegerVector idx,int len){
  int nbin=((len-1)>>5)+1;
  int i;
  IntegerVector output(nbin);
  for(i=0;i<idx.size();i++){
    int bin_x=idx(i)>>5;
    int bin_y=idx(i)&31;
    output(bin_x)=output(bin_x)|(1<<bin_y);
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector locusAND(IntegerVector ve1,IntegerVector ve2){
  int i;
  IntegerVector output(ve1.size());
  for(i=0;i<ve1.size();i++){
    output(i)=ve1(i)&ve2(i);
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector locusXOR(IntegerVector ve1,IntegerVector ve2){
  int i;
  IntegerVector output(ve1.size());
  for(i=0;i<ve1.size();i++){
    output(i)=ve1(i)^ve2(i);
  }
  return output;
}

// [[Rcpp::export]]
double locusSUM(IntegerVector ve,NumericVector weight){
  double output=0;
  for(int x=0;x<ve.size();x++){
    int mask=ve(x);
    for(int y=0;y<32&&x*32+y<weight.size();y++){
      if(((mask>>y)&1)==1){
        output+=weight(x*32+y);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
int locusSUMint(IntegerVector ve){
  int output=0;
  for(int x=0;x<ve.size();x++){
    int mask=ve(x);
    for(int y=0;y<32;y++){
      if(((mask>>y)&1)==1){
        output++;
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix reOrderEdge(IntegerMatrix edge,int root){
  IntegerMatrix output(edge.nrow(),edge.ncol());
  IntegerVector ances_stack(edge.nrow()+1);
  ances_stack(0)=root;
  int stack_pointer=0;
  bool isfind;
  int i=0,j;
  while(stack_pointer>(-1)){
    isfind=false;
    for(j=0;j<edge.nrow();j++){
      if(edge(j,0)==ances_stack(stack_pointer)){
        output.row(i)=edge.row(j);
        i++;
        stack_pointer++;
        ances_stack(stack_pointer)=edge(j,1);
        edge(j,0)=-1;
        isfind=true;
        break;
      }
    }
    if(!isfind){
      stack_pointer--;
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector charAncestral(IntegerMatrix mat){
  int i,j;
  IntegerVector output(mat.ncol());
  for(i=0;i<mat.ncol();i++){
    output(i)=(-1);
    for(j=0;j<mat.nrow();j++){
      if(mat(j,i)>=0&&output(i)==-1){
        output(i)=mat(j,i);
      }else if(mat(j,i)>=0&&mat(j,i)!=output(i)){
        output(i)=0;
        break;
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector bitAncestral(IntegerMatrix mat){
  int i,j;
  IntegerVector output(mat.ncol());
  for(j=0;j<mat.ncol();j++){
    output(j)=mat(0,j);
  }
  for(i=1;i<mat.nrow();i++){
    for(j=0;j<mat.ncol();j++){
      output(j)=output(j)&mat(i,j);
    }
  }
  return output;
}

double countApo(IntegerVector ve1,IntegerVector ve2,NumericMatrix weight){
  double nshared=0;
  int i;
  for(i=0;i<ve1.size();i++){
    if(ve1(i)==ve2(i)&&ve1(i)>0){
      nshared+=weight(ve1(i)-1,i);
    }
  }
  return nshared;
}

float countApoFloat(IntegerVector ve1,IntegerVector ve2,NumericMatrix weight){
  float nshared=0;
  int i;
  for(i=0;i<ve1.size();i++){
    if(ve1(i)==ve2(i)&&ve1(i)>0){
      nshared+=weight(ve1(i)-1,i);
    }
  }
  return nshared;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdge(IntegerMatrix mat,NumericMatrix weight,bool verbose=false){
  IntegerMatrix edge(mat.nrow()*2-2,3+mat.ncol());
  IntegerVector node2cache(mat.nrow()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.nrow());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  float** cacheShared=(float**)(malloc(mat.nrow()*sizeof(float*)));
  for(int i=0;i<mat.nrow();i++){
    cacheShared[i]=(float*)(malloc(mat.nrow()*sizeof(float)));
  }
  float rowMax[mat.nrow()*2-1];
  for(int i=0;i<mat.nrow()*2-1;i++){
    rowMax[i]=-1;
  }
  IntegerVector rowMaxIdx(mat.nrow()*2-1);
  
  for(int i=0;i<mat.nrow();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.row(i)=mat.row(i);
  }
  for(int i=0;i<mat.nrow()-1;i++){
    for(int j=i+1;j<mat.nrow();j++){
      float shared=countApoFloat(mat.row(i),mat.row(j),weight);
      cacheShared[i][j]=shared;
      if(rowMax[i]<shared){
        rowMax[i]=shared;
        rowMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.nrow();
  while(nodeCount<mat.nrow()*2-1){
    if(verbose&&((nodeCount-mat.nrow())&1023)==1023){
      Rprintf("=");
    }
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(rowMaxIdx(i))==-1){
          rowMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(rowMax[i]<shared){
                rowMax[i]=shared;
                rowMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<rowMax[i]){
          maxShared=rowMax[i];
          joinIdx1=i;
          joinIdx2=rowMaxIdx(i);
        }
      }
    }
    
    edge((nodeCount-mat.nrow())*2,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2,1)=joinIdx1;
    edge((nodeCount-mat.nrow())*2+1,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2+1,1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerMatrix tmpMat(2,mat.ncol());
    tmpMat.row(0)=cacheMat.row(cacheIdx1);
    tmpMat.row(1)=cacheMat.row(cacheIdx2);
    IntegerVector stateAns=charAncestral(tmpMat);
    
    for(int k=0;k<mat.ncol();k++){
      edge((nodeCount-mat.nrow())*2,3+k)=cacheMat(cacheIdx1,k);
      edge((nodeCount-mat.nrow())*2+1,3+k)=cacheMat(cacheIdx2,k);
      if(stateAns(k)==0&&cacheMat(cacheIdx1,k)>0){
        edge((nodeCount-mat.nrow())*2,2)++;
      }
      if(stateAns(k)==0&&cacheMat(cacheIdx2,k)>0){
        edge((nodeCount-mat.nrow())*2+1,2)++;
      }
    }
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.row(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=countApoFloat(cacheMat.row(tmp),cacheMat.row(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(rowMax[i]<shared){
          rowMax[i]=shared;
          rowMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

float countApoBit(IntegerVector ve1,IntegerVector ve2,NumericVector weight){
  float output=0;
  for(int x=0;x<ve1.size();x++){
    int mask=ve1(x)&ve2(x);
    for(int y=0;y<32&&x*32+y<weight.size();y++){
      if(((mask>>y)&1)==1){
        output+=weight(x*32+y);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdgeBit(IntegerMatrix mat,NumericVector weight,bool verbose=false){
  IntegerMatrix edge(mat.nrow()*2-2,3+mat.ncol());
  IntegerVector node2cache(mat.nrow()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.nrow());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  float** cacheShared=(float**)(malloc(mat.nrow()*sizeof(float*)));
  for(int i=0;i<mat.nrow();i++){
    cacheShared[i]=(float*)(malloc(mat.nrow()*sizeof(float)));
  }
  float rowMax[mat.nrow()*2-1];
  for(int i=0;i<mat.nrow()*2-1;i++){
    rowMax[i]=-1;
  }
  IntegerVector rowMaxIdx(mat.nrow()*2-1);
  
  for(int i=0;i<mat.nrow();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.row(i)=mat.row(i);
  }
  for(int i=0;i<mat.nrow()-1;i++){
    for(int j=i+1;j<mat.nrow();j++){
      float shared=countApoBit(mat.row(i),mat.row(j),weight=weight);
      cacheShared[i][j]=shared;
      if(rowMax[i]<shared){
        rowMax[i]=shared;
        rowMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.nrow();
  while(nodeCount<mat.nrow()*2-1){
    if(verbose&&((nodeCount-mat.nrow())&1023)==1023){
      Rprintf("=");
    }
    
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(rowMaxIdx(i))==-1){
          rowMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(rowMax[i]<shared){
                rowMax[i]=shared;
                rowMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<rowMax[i]){
          maxShared=rowMax[i];
          joinIdx1=i;
          joinIdx2=rowMaxIdx(i);
        }
      }
    }
    
    edge((nodeCount-mat.nrow())*2,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2,1)=joinIdx1;
    edge((nodeCount-mat.nrow())*2+1,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2+1,1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerMatrix tmpMat(2,mat.ncol());
    tmpMat.row(0)=cacheMat.row(cacheIdx1);
    tmpMat.row(1)=cacheMat.row(cacheIdx2);
    IntegerVector stateAns=bitAncestral(tmpMat);
    
    for(int k=0;k<mat.ncol();k++){
      edge((nodeCount-mat.nrow())*2,3+k)=cacheMat(cacheIdx1,k);
      edge((nodeCount-mat.nrow())*2+1,3+k)=cacheMat(cacheIdx2,k);
    }
    
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.row(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=countApoBit(cacheMat.row(tmp),cacheMat.row(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(rowMax[i]<shared){
          rowMax[i]=shared;
          rowMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

// [[Rcpp::export]]
IntegerMatrix collapseBranch(IntegerMatrix edge,int ntip){
  int i,j;
  for(i=0;i<edge.nrow();i++){
    if(edge(i,2)==0&&edge(i,1)>=ntip){
      for(j=0;j<edge.nrow();j++){
        if(edge(j,0)==edge(i,1)){
          edge(j,0)=edge(i,0);
        }
        edge(i,2)=-1;
      }
    }
  }
  return edge;
}

// [[Rcpp::export]]
IntegerVector nodeKeep(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
  int i,j,k=0;
  int Nnode=node_mat.nrow();
  IntegerVector state(Ntip+Nnode);
  IntegerVector output(Ntip+Nnode);
  state.fill(-1);
  for(i=0;i<Ntip;i++){
    state(i)=1;
  }
  for(i=edge.nrow()-1;i>=0;i--){
    if(state(edge(i,1)-1)==0){
      state(edge(i,0)-1)=0;
    }else if(state(edge(i,1)-1)==1){
      if(state(edge(i,0)-1)==-1){
        state(edge(i,0)-1)=1;
        for(j=0;j<outmat.nrow();j++){
          if(mutationNested(outmat.row(j),node_mat.row(edge(i,0)-1-Ntip))<keep_cutoff){
            state(edge(i,0)-1)=0;
            break;
          }
        }
      }
      if(state(edge(i,0)-1)==0){
        output(k)=edge(i,1);
        k++;
      }
    }else if(state(edge(i,1)-1)==-1){
      Rprintf("Error at row %d\n",i+1);
      return state;
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector nodeKeepBit(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
  int i,j,k=0;
  int Nnode=node_mat.nrow();
  IntegerVector state(Ntip+Nnode);
  IntegerVector output(Ntip+Nnode);
  state.fill(-1);
  for(i=0;i<Ntip;i++){
    state(i)=1;
  }
  for(i=edge.nrow()-1;i>=0;i--){
    if(state(edge(i,1)-1)==0){
      state(edge(i,0)-1)=0;
    }else if(state(edge(i,1)-1)==1){
      if(state(edge(i,0)-1)==-1){
        state(edge(i,0)-1)=1;
        for(j=0;j<outmat.nrow();j++){
          if(mutationNestedBit(outmat.row(j),node_mat.row(edge(i,0)-1-Ntip))<keep_cutoff){
            state(edge(i,0)-1)=0;
            break;
          }
        }
      }
      if(state(edge(i,0)-1)==0){
        output(k)=edge(i,1);
        k++;
      }
    }else if(state(edge(i,1)-1)==-1){
      Rprintf("Error at row %d\n",i+1);
      return state;
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector imputeVector(IntegerVector ref,IntegerVector query){
  int i;
  IntegerVector output(ref.size());
  for(i=0;i<ref.size();i++){
    if(query(i)==-1){
      output(i)=ref(i);
    }else{
      output(i)=query(i);
    }
  }
  return output;
}

// [[Rcpp::export]]
int markDup(int row_idx,IntegerMatrix mat){
  int output=row_idx;
  int j,x,y;
  unsigned int k;
  std::vector<int> idx;
  std::vector<int> xidx;
  std::vector<int> yidx;
  for(x=0;x<mat.ncol();x++){
    for(y=0;y<32;y++){
      if(((mat(row_idx,x)>>y)&1)==1){
        idx.push_back(x*32+y);
        xidx.push_back(x);
        yidx.push_back(y);
      }
    }
  }
  unsigned nshared;
  for(j=0;j<mat.nrow();j++){
    nshared=0;
    for(k=0;k<idx.size();k++){
      if(((mat(j,xidx[k])>>yidx[k])&1)==1){
        nshared++;
      }else{
        break;
      }
    }
    if(nshared==idx.size()&&output!=j){
      output=j;
      break;
    }
  }
  if(output<=row_idx){
    return output;
  }else{
    unsigned count_output=0;
    for(x=0;x<mat.ncol();x++){
      for(y=0;y<32;y++){
        if(((mat(output,x)>>y)&1)==1){
          count_output++;
        }
      }
    }
    if(count_output>idx.size()){
      return output;
    }else{
      return row_idx;
    }
  }
}

// [[Rcpp::export]]
int markDupBit(int row_idx,IntegerMatrix mat){
  int output=row_idx;
  int j=-1;
  bool dup=false;
  
  do{
    j++;
    output=j;
    dup=true;
    for(int k=0;k<mat.ncol();k++){
      if(mat(row_idx,k)!=mat(j,k)){
        dup=false;
        break;
      }
    }
  } while (!dup);
  return output;
}

// [[Rcpp::export]]
IntegerMatrix edgeEvo(IntegerMatrix edge,IntegerMatrix binary,int root){
  IntegerMatrix output(binary.nrow(),binary.ncol());
  int i,j;
  for(i=0;i<edge.nrow();i++){
    for(j=0;j<binary.ncol();j++){
      output(edge(i,1)-1,j)=((~binary(edge(i,0)-1,j))&binary(edge(i,1)-1,j));
    }
  }
  output.row(root)=binary.row(root);
  return output;
}

// [[Rcpp::export]]
void printB(int x){
  int i;
  for(i=31;i>=0;i--){
    Rprintf("%d",(x>>i)&1);
  }
  Rprintf("\n");
}

// [[Rcpp::export]]
IntegerVector edge2vector(IntegerMatrix edge,int maxIdx){
  int i;
  IntegerVector route(maxIdx);
  route.fill(-1);
  for(i=0;i<edge.nrow();i++){
    route(edge(i,1)-1)=edge(i,0)-1;
  }
  return route;
}

template <typename T>

IntegerVector idxSort(T* arr,int len){
  IntegerVector idx(len);
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),[&](int i,int j){return arr[i]>arr[j];});
  return idx;
}

// [[Rcpp::export]]
IntegerMatrix findNeighborWeighted(IntegerMatrix mat,int maxChar,NumericVector weight,int n=10,double p=0.5,double radius=0,bool verbose=false){
  int i,j,k,x,y;
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerMatrix mat_neighbor(mat.nrow(),nbin);
  if(verbose){
    Rprintf("Finding neighbors of matrix with %d rows.\n",mat.nrow());
  }
  for(i=0;i<mat.nrow();i++){
    if(verbose&&(i&1023)==0){
      Rprintf("=");
    }
    int idx[maxChar];
    int xidx[maxChar];
    int yidx[maxChar];
    int nchar=0;
    for(x=0;x<mat.ncol();x++){
      for(y=0;y<32;y++){
        if(((mat(i,x)>>y)&1)==1){
          idx[nchar]=x*32+y;
          xidx[nchar]=x;
          yidx[nchar]=y;
          nchar++;
        }
      }
    }
    if(radius!=0){
      double nshared=0;
      for(k=0;k<nchar;k++){
        nshared+=weight(idx[k]);
      }
      if(nshared>=radius){
        int bin_x=i>>5;
        int bin_y=i&31;
        mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
      }
      for(j=i+1;j<mat.nrow();j++){
        nshared=0;
        for(k=0;k<nchar;k++){
          if(((mat(j,xidx[k])>>yidx[k])&1)==1){
            nshared+=weight(idx[k]);
          }
        }
        if(nshared>=radius){
          int bin_x=j>>5;
          int bin_y=j&31;
          mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
          bin_x=i>>5;
          bin_y=i&31;
          mat_neighbor(j,bin_x)=mat_neighbor(j,bin_x)|(1<<bin_y);
        }
      }
    }else{
      double dist[mat.nrow()];
      int distIdx[mat.nrow()];
      int nDistIdx=0;
      double cutoff=nchar*p;
      for(j=0;j<mat.nrow();j++){
        double nshared=0;
        for(k=0;k<nchar;k++){
          if(((mat(j,xidx[k])>>yidx[k])&1)==1){
            nshared+=weight(idx[k]);
          }
        }
        if(nshared>=cutoff){
          dist[nDistIdx]=nshared;
          distIdx[nDistIdx]=j;
          nDistIdx++;
        }
      }
      IntegerVector sortedIdx=idxSort(dist,nDistIdx);
      if(sortedIdx.size()>n){
        cutoff=dist[sortedIdx(n-1)];
      }
      for(j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
        int bin_x=distIdx[sortedIdx(j)]>>5;
        int bin_y=distIdx[sortedIdx(j)]&31;
        mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
      }
    }
  }
  if(verbose){
    Rprintf("\n");
  }
  return mat_neighbor;
}

// [[Rcpp::export]]
IntegerMatrix findNeighborUnweighted(IntegerMatrix mat,int maxChar,int n=10,double p=0.5,int radius=0,bool verbose=false){
  int i,j,k,x,y;
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerMatrix mat_neighbor(mat.nrow(),nbin);
  if(verbose){
    Rprintf("Finding neighbors of matrix with %d rows.\n",mat.nrow());
  }
  for(i=0;i<mat.nrow();i++){
    if(verbose&&(i&1023)==0){
      Rprintf("=");
    }
    int xidx[maxChar];
    int yidx[maxChar];
    int nchar=0;
    for(x=0;x<mat.ncol();x++){
      for(y=0;y<32;y++){
        if(((mat(i,x)>>y)&1)==1){
          xidx[nchar]=x;
          yidx[nchar]=y;
          nchar++;
        }
      }
    }
    if(radius!=0){
      int nshared=0;
      for(k=0;k<nchar;k++){
        nshared++;
      }
      if(nshared>=radius){
        int bin_x=i>>5;
        int bin_y=i&31;
        mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
      }
      for(j=i+1;j<mat.nrow();j++){
        nshared=0;
        for(k=0;k<nchar;k++){
          if(((mat(j,xidx[k])>>yidx[k])&1)==1){
            nshared++;
          }
        }
        if(nshared>=radius){
          int bin_x=j>>5;
          int bin_y=j&31;
          mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
          bin_x=i>>5;
          bin_y=i&31;
          mat_neighbor(j,bin_x)=mat_neighbor(j,bin_x)|(1<<bin_y);
        }
      }
    }else{
      int dist[mat.nrow()];
      int distIdx[mat.nrow()];
      int nDistIdx=0;
      int cutoff=static_cast<int>(nchar*p);
      for(j=0;j<mat.nrow();j++){
        int nshared=0;
        for(k=0;k<nchar;k++){
          if(((mat(j,xidx[k])>>yidx[k])&1)==1){
            nshared++;
          }
        }
        if(nshared>=cutoff){
          dist[nDistIdx]=nshared;
          distIdx[nDistIdx]=j;
          nDistIdx++;
        }
      }
      IntegerVector sortedIdx=idxSort(dist,nDistIdx);
      if(sortedIdx.size()>n){
        cutoff=dist[sortedIdx(n-1)];
      }
      for(j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
        int bin_x=distIdx[sortedIdx(j)]>>5;
        int bin_y=distIdx[sortedIdx(j)]&31;
        mat_neighbor(i,bin_x)=mat_neighbor(i,bin_x)|(1<<bin_y);
      }
    }
  }
  if(verbose){
    Rprintf("\n");
  }
  return mat_neighbor;
}

// [[Rcpp::export]]
IntegerVector findNeighborWeightedParallel(int i,IntegerMatrix mat,int maxChar,NumericVector weight,int n=10,double p=0.5,double radius=0){
  int j,k,x,y;
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerVector mat_neighbor(nbin);
  int idx[maxChar];
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  for(x=0;x<mat.ncol();x++){
    for(y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        idx[nchar]=x*32+y;
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
      }
    }
  }
  if(radius!=0){
    for(j=0;j<mat.nrow();j++){
      double nshared=0;
      for(k=0;k<nchar;k++){
        if(((mat(j,xidx[k])>>yidx[k])&1)==1){
          nshared+=weight(idx[k]);
        }
      }
      if(nshared>=radius){
        int bin_x=j>>5;
        int bin_y=j&31;
        mat_neighbor(bin_x)|=(1<<bin_y);
      }
    }
  }else{
    double dist[mat.nrow()];
    int distIdx[mat.nrow()];
    int nDistIdx=0;
    double cutoff=nchar*p;
    for(j=0;j<mat.nrow();j++){
      double nshared=0;
      for(k=0;k<nchar;k++){
        if(((mat(j,xidx[k])>>yidx[k])&1)==1){
          nshared+=weight(idx[k]);
        }
      }
      if(nshared>=cutoff){
        dist[nDistIdx]=nshared;
        distIdx[nDistIdx]=j;
        nDistIdx++;
      }
    }
    IntegerVector sortedIdx=idxSort(dist,nDistIdx);
    if(sortedIdx.size()>n){
      cutoff=dist[sortedIdx(n-1)];
    }
    for(j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
      int bin_x=distIdx[sortedIdx(j)]>>5;
      int bin_y=distIdx[sortedIdx(j)]&31;
      mat_neighbor(bin_x)=mat_neighbor(bin_x)|(1<<bin_y);
    }
  }
  return mat_neighbor;
}

// [[Rcpp::export]]
IntegerVector findNeighborUnweightedParallel(int i,IntegerMatrix mat,int maxChar,int n=10,double p=0.5,int radius=0){
  int j,k,x,y;
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerVector mat_neighbor(nbin);
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  for(x=0;x<mat.ncol();x++){
    for(y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
      }
    }
  }
  if(radius!=0){
    for(j=0;j<mat.nrow();j++){
      int nshared=0;
      for(k=0;k<nchar;k++){
        if(((mat(j,xidx[k])>>yidx[k])&1)==1){
          nshared++;
        }
      }
      if(nshared>=radius){
        int bin_x=j>>5;
        int bin_y=j&31;
        mat_neighbor(bin_x)|=(1<<bin_y);
      }
    }
  }else{
    int dist[mat.nrow()];
    int distIdx[mat.nrow()];
    int nDistIdx=0;
    int cutoff=static_cast<int>(nchar*p);
    for(j=0;j<mat.nrow();j++){
      int nshared=0;
      for(k=0;k<nchar;k++){
        if(((mat(j,xidx[k])>>yidx[k])&1)==1){
          nshared++;
        }
      }
      if(nshared>=cutoff){
        dist[nDistIdx]=nshared;
        distIdx[nDistIdx]=j;
        nDistIdx++;
      }
    }
    IntegerVector sortedIdx=idxSort(dist,nDistIdx);
    if(sortedIdx.size()>n){
      cutoff=dist[sortedIdx(n-1)];
    }
    for(j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
      int bin_x=distIdx[sortedIdx(j)]>>5;
      int bin_y=distIdx[sortedIdx(j)]&31;
      mat_neighbor(bin_x)=mat_neighbor(bin_x)|(1<<bin_y);
    }
  }
  return mat_neighbor;
}

// [[Rcpp::export]]
IntegerVector tipMergeNeighbor(IntegerMatrix mat,double p=0.8,bool sqrt_nNeighbor=false,
                                int min_shared=0,bool verbose=false){
  int i,j,k,x,y;
  IntegerVector nNeighbor(mat.nrow());
  for(i=0;i<mat.nrow();i++){
    for(x=0;x<mat.ncol();x++){
      for(y=0;y<32;y++){
        if(((mat(i,x)>>y)&1)==1){
          nNeighbor(i)++;
        }
      }
    }
  }
  int cmpCount=0;
  IntegerVector member(mat.nrow());
  member.fill(-1);
  IntegerVector transition(mat.nrow());
  transition.fill(-1);
  int nCategory=0;
  if(verbose){
    Rprintf("Clustering %d tips.\n",mat.nrow());
  }
  for(i=0;i<mat.nrow();i++){
    if(verbose&&(i&1023)==1023){
      Rprintf("=");
    }
    int root1;
    if(member(i)==-1){
      root1=nCategory;
      member(i)=nCategory;
      transition(nCategory)=nCategory;
      nCategory++;
    }else{
      root1=transitionRoot(transition,member(i));
    }
    int nCount=0;
    int neighbor_idx[nNeighbor(i)];
    int neighbor_xidx[nNeighbor(i)];
    int neighbor_yidx[nNeighbor(i)];
    for(x=0;x<mat.ncol();x++){
      for(y=0;y<32&&nCount<nNeighbor(i);y++){
        if(((mat(i,x)>>y)&1)==1){
          neighbor_idx[nCount]=x*32+y;
          neighbor_xidx[nCount]=x;
          neighbor_yidx[nCount]=y;
          nCount++;
        }
      }
    }
    for(j=0;j<nCount;j++){
      if(i>=neighbor_idx[j]){
        continue;
      }
      int current_idx=neighbor_idx[j];
      int root2=transitionRoot(transition,member(current_idx));
      if(root1!=root2){
        cmpCount++;
        int nshared=0;
        for(k=0;k<nCount;k++){
          if(((mat(current_idx,neighbor_xidx[k])>>neighbor_yidx[k])&1)==1){
            nshared++;
          }
        }
        if(nshared<min_shared){
          continue;
        }
        double expected;
        if(sqrt_nNeighbor){
          expected=sqrt(nCount*nNeighbor(current_idx))*p;
        }else{
          IntegerVector cmp_temp={nCount,nNeighbor(current_idx)};
          expected=min(cmp_temp)*p;
        }
        if(nshared>expected){
          if(root2==-1){
            member(current_idx)=root1;
          }else{
            transition(root1)=root2;
            root1=root2;
          }
        }
      }
    }
  }
  for(i=0;i<member.size();i++){
    member(i)=transitionRoot(transition,member(i))+1;
  }
  if(verbose){
    Rprintf("\n");
  }
  return member;
}

// [[Rcpp::export]]
IntegerVector tipMergeNeighborBig(SEXP pNeighborMat,IntegerVector nNeighbor,
                                  double p=0.8,bool sqrt_nNeighbor=false,
                                  int min_shared=0,bool verbose=false){
  XPtr<BigMatrix> xpMat(pNeighborMat);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int i,j,k,x,y;
  int root1,root2,nCount,mask,current_idx,nshared;
  int xidx,xidx2,yidx;
  double expected;
  int n_cols=xpMat->ncol();
  int n_rows=xpMat->nrow();
  
  IntegerVector member(n_cols);
  member.fill(-1);
  IntegerVector transition(n_cols);
  transition.fill(-1);
  int nCategory=0;
  if(verbose){
    Rprintf("Clustering %d tips.\n",n_cols);
  }
  int* neighbor_idx=(int*)(malloc(n_cols*sizeof(int)));
  for(i=0;i<n_cols;i++){
    if(verbose&&(i&1023)==1023){
      Rprintf("=");
    }
    if(member(i)==-1){
      root1=nCategory;
      member(i)=nCategory;
      transition(nCategory)=nCategory;
      nCategory++;
    }else{
      root1=transitionRoot(transition,member(i));
    }
    nCount=0;
    for(x=0;x<n_rows;x++){
      mask=mat[i][x];
      for(y=0;y<32&&nCount<nNeighbor(i);y++){
        if(((mask>>y)&1)==1){
          neighbor_idx[nCount]=x*32+y;
          nCount++;
        }
      }
    }
    for(j=0;j<nCount;j++){
      current_idx=neighbor_idx[j];
      if(i>=current_idx){
        continue;
      }
      root2=transitionRoot(transition,member(current_idx));
      if(root1!=root2){
        nshared=0;
        xidx=-1;
        mask=0;
        for(k=0;k<nCount;k++){
          yidx=neighbor_idx[k]&31;
          xidx2=neighbor_idx[k]>>5;
          if(xidx!=xidx2){
            xidx=xidx2;
            mask=mat[current_idx][xidx];
          }
          if(((mask>>yidx)&1)==1){
            nshared++;
          }
        }
        if(nshared<min_shared){
          continue;
        }
        if(sqrt_nNeighbor){
          expected=sqrt(nCount*nNeighbor(current_idx))*p;
        }else{
          if(nCount<nNeighbor(current_idx)){
            expected=nCount*p;
          }else{
            expected=nNeighbor(current_idx)*p;
          }
        }
        if(nshared>expected){
          if(root2==-1){
            member(current_idx)=root1;
          }else{
            transition(root1)=root2;
            root1=root2;
          }
        }
      }
    }
  }
  for(i=0;i<member.size();i++){
    member(i)=transitionRoot(transition,member(i))+1;
  }
  if(verbose){
    Rprintf("\n");
  }
  return member;
}

// [[Rcpp::export]]
IntegerMatrix nodeDes(IntegerMatrix edge,int nnode,int ntip){
  int nbin=((ntip-1)>>5)+1;
  IntegerMatrix output(nbin,nnode);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1-ntip;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      output(x,ansIdx)|=(1<<y);
    }else{
      for(int j=0;j<output.nrow();j++){
        output(j,ansIdx)|=output(j,desIdx-ntip);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
void nodeDesBig(IntegerMatrix edge,int nnode,int ntip,std::string file_path,int cacheSize=8192){
  int nbin=((ntip-1)>>5)+1;
  int** cacheMat=(int**)(malloc(cacheSize*sizeof(int*)));
  for(int i=0;i<cacheSize;i++){
    cacheMat[i]=(int*)(malloc(nbin*sizeof(int)));
  }
  for(int i=0;i<cacheSize;i++){
    for(int j=0;j<nbin;j++){
      cacheMat[i][j]=0;
    }
  }
  LogicalVector cacheFree(cacheSize);
  cacheFree.fill(true);
  int nUse=0;
  IntegerVector node2cache(nnode);
  node2cache.fill(-1);
  int cacheIdx=0;
  int dim[2];
  dim[0]=nnode-1;dim[1]=nbin;

  FILE *file;
  file = fopen(file_path.c_str(), "wb");
  fwrite(dim,sizeof(int),2,file);
  
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      int ansCacheIdx=node2cache(ansIdx-ntip);
      if(ansCacheIdx==-1){
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return;
        }
        while(!cacheFree(cacheIdx)){
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        ansCacheIdx=cacheIdx;
        node2cache(ansIdx-ntip)=cacheIdx;
        cacheFree(cacheIdx)=false;
        nUse++;
      }
      cacheMat[ansCacheIdx][x]|=(1<<y);
    }else{
      if(node2cache(ansIdx-ntip)==-1){
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return;
        }
        while(!cacheFree(cacheIdx)){
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        node2cache(ansIdx-ntip)=cacheIdx;
        cacheFree(cacheIdx)=false;
        nUse++;
      }
      int desCacheIdx=node2cache(desIdx-ntip);
      int ansCacheIdx=node2cache(ansIdx-ntip);
      for(int x=0;x<nbin;x++){
        cacheMat[ansCacheIdx][x]|=cacheMat[desCacheIdx][x];
      }
      fwrite(cacheMat[desCacheIdx],sizeof(int),nbin,file);
      cacheFree(desCacheIdx)=true;
      nUse--;
      for(int x=0;x<nbin;x++){
        cacheMat[desCacheIdx][x]=0;
      }
    }
  }
  fclose(file);
  return;
}

// [[Rcpp::export]]
IntegerVector findNeighborRangeUnweighted(int i,IntegerMatrix mat,int maxChar,IntegerMatrix range, int offset, int n=10,double p=0.5,int radius=0){
  IntegerVector output(range.ncol());
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  for(int x=0;x<mat.ncol();x++){
    for(int y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
      }
    }
  }
  if(radius!=0){
    for(int x=0;x<range.ncol();x++){
      for(int y=0;x*32+y<mat.nrow()&&y<32;y++){
        if(((range(i-offset,x)>>y)&1)==1){
          int jIdx=x*32+y;
          int nshared=0;
          for(int k=0;k<nchar;k++){
            if(((mat(jIdx,xidx[k])>>yidx[k])&1)==1){
              nshared++;
            }
          }
          if(nshared>=radius){
            output(x)|=(1<<y);
          }
        }
      }
    }
  }else{
    int dist[mat.nrow()];
    int distIdx[mat.nrow()];
    int nDistIdx=0;
    int cutoff=static_cast<int>(nchar*p);
    for(int x=0;x<range.ncol();x++){
      for(int y=0;x*32+y<mat.nrow()&&y<32;y++){
        if(((range(i-offset,x)>>y)&1)==1){
          int jIdx=x*32+y;
          int nshared=0;
          for(int k=0;k<nchar;k++){
            if(((mat(jIdx,xidx[k])>>yidx[k])&1)==1){
              nshared++;
            }
          }
          if(nshared>=cutoff){
            dist[nDistIdx]=nshared;
            distIdx[nDistIdx]=jIdx;
            nDistIdx++;
          }
        }
      }
    }
    IntegerVector sortedIdx=idxSort(dist,nDistIdx);
    if(sortedIdx.size()>n){
      cutoff=dist[sortedIdx(n-1)];
    }
    for(int j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
      int tmp=distIdx[sortedIdx(j)];
      output(tmp>>5)|=(1<<(tmp&31));
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector findNeighborRangeWeighted(int i,IntegerMatrix mat,int maxChar,IntegerMatrix range, int offset,NumericVector weight,int n=10,double p=0.5,double radius=0){
  IntegerVector output(range.ncol());
  int idx[maxChar];
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  double radiusSelf=0;
  for(int x=0;x<mat.ncol();x++){
    for(int y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        idx[nchar]=x*32+y;
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
        radiusSelf+=weight(x*32+y);
      }
    }
  }
  if(radius!=0){
    for(int x=0;x<range.ncol();x++){
      for(int y=0;x*32+y<mat.nrow()&&y<32;y++){
        if(((range(i-offset,x)>>y)&1)==1){
          int jIdx=x*32+y;
          double nshared=0;
          for(int k=0;k<nchar;k++){
            if(((mat(jIdx,xidx[k])>>yidx[k])&1)==1){
              nshared+=weight(idx[k]);
            }
          }
          if(nshared>=radius){
            output(x)|=(1<<y);
          }
        }
      }
    }
  }else{
    double dist[mat.nrow()];
    int distIdx[mat.nrow()];
    int nDistIdx=0;
    double cutoff=radiusSelf*p;
    for(int x=0;x<range.ncol();x++){
      for(int y=0;x*32+y<mat.nrow()&&y<32;y++){
        if(((range(i-offset,x)>>y)&1)==1){
          int jIdx=x*32+y;
          double nshared=0;
          for(int k=0;k<nchar;k++){
            if(((mat(jIdx,xidx[k])>>yidx[k])&1)==1){
              nshared+=weight(idx[k]);
            }
          }
          if(nshared>=cutoff){
            dist[nDistIdx]=nshared;
            distIdx[nDistIdx]=jIdx;
            nDistIdx++;
          }
        }
      }
    }
    IntegerVector sortedIdx=idxSort(dist,nDistIdx);
    if(sortedIdx.size()>n){
      cutoff=dist[sortedIdx(n-1)];
    }
    for(int j=0;j<sortedIdx.size()&&dist[sortedIdx(j)]>=cutoff;j++){
      int tmp=distIdx[sortedIdx(j)];
      output(tmp>>5)|=(1<<(tmp&31));
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector imputeMat(int i,IntegerMatrix mat, IntegerMatrix neighbor_mat,int offset) {
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  int neighbor_ncols = neighbor_mat.ncol();
  IntegerVector output=mat.row(i);
  
  std::vector<int> mask_loci;
  for (int j = 0; j < ncols; ++j) {
    if (mat(i, j) == -1) {
      mask_loci.push_back(j);
    }
  }
  
  if (!mask_loci.empty()) {
    std::vector<int> neighbor;
    for (int x = 0; x < neighbor_ncols; ++x) {
      int mask = neighbor_mat(i-offset, x);
      for (int y = 0; y < 32; ++y) {
        if ((mask >> y) & 1) { // Check if the y-th bit is set
          int neighbor_idx = x * 32 + y;
          if (neighbor_idx < nrows) { // Ensure neighbor index is valid
            neighbor.push_back(neighbor_idx);
          }
        }
      }
    }
    
    for (size_t j_idx = 0; j_idx < mask_loci.size(); ++j_idx) {
      int j = mask_loci[j_idx];
      int state_output = -1;
      for (size_t k = 0; k < neighbor.size(); ++k) {
        int n_row = neighbor[k];
        int state = mat(n_row, j);
        if (state >= 0) {
          if (state_output == -1) {
            state_output = state;
          } else if (state != state_output) {
            state_output = 0;
            break;
          }
        }
      }
      if (state_output > 0) {
        output(j) = state_output;
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdgeBig(IntegerMatrix mat,NumericMatrix weight,SEXP pShared,bool verbose=false){
  XPtr<BigMatrix> xpShared(pShared);
  MatrixAccessor<double> cacheShared=MatrixAccessor<double>(*xpShared);
  IntegerMatrix edge(mat.nrow()*2-2,3+mat.ncol());
  IntegerVector node2cache(mat.nrow()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.nrow());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  double rowMax[mat.nrow()*2-1];
  for(int i=0;i<mat.nrow()*2-1;i++){
    rowMax[i]=-1;
  }
  IntegerVector rowMaxIdx(mat.nrow()*2-1);
  
  for(int i=0;i<mat.nrow();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.row(i)=mat.row(i);
  }
  for(int i=0;i<mat.nrow()-1;i++){
    for(int j=i+1;j<mat.nrow();j++){
      double shared=countApo(mat.row(i),mat.row(j),weight);
      cacheShared[i][j]=shared;
      if(rowMax[i]<shared){
        rowMax[i]=shared;
        rowMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.nrow();
  while(nodeCount<mat.nrow()*2-1){
    if(verbose&&((nodeCount-mat.nrow())&1023)==1023){
      Rprintf("=");
    }
    double maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(rowMaxIdx(i))==-1){
          rowMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              double shared=cacheShared[cacheIdxi][cacheIdxj];
              if(rowMax[i]<shared){
                rowMax[i]=shared;
                rowMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<rowMax[i]){
          maxShared=rowMax[i];
          joinIdx1=i;
          joinIdx2=rowMaxIdx(i);
        }
      }
    }
    
    edge((nodeCount-mat.nrow())*2,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2,1)=joinIdx1;
    edge((nodeCount-mat.nrow())*2+1,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2+1,1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerMatrix tmpMat(2,mat.ncol());
    tmpMat.row(0)=cacheMat.row(cacheIdx1);
    tmpMat.row(1)=cacheMat.row(cacheIdx2);
    IntegerVector stateAns=charAncestral(tmpMat);
    
    for(int k=0;k<mat.ncol();k++){
      edge((nodeCount-mat.nrow())*2,3+k)=cacheMat(cacheIdx1,k);
      edge((nodeCount-mat.nrow())*2+1,3+k)=cacheMat(cacheIdx2,k);
    }
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.row(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        double shared=countApo(cacheMat.row(tmp),cacheMat.row(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(rowMax[i]<shared){
          rowMax[i]=shared;
          rowMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdgeBitBig(IntegerMatrix mat,NumericVector weight,SEXP pShared,bool verbose=false){
  XPtr<BigMatrix> xpShared(pShared);
  MatrixAccessor<double> cacheShared=MatrixAccessor<double>(*xpShared);
  IntegerMatrix edge(mat.nrow()*2-2,3+mat.ncol());
  IntegerVector node2cache(mat.nrow()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.nrow());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  double rowMax[mat.nrow()*2-1];
  for(int i=0;i<mat.nrow()*2-1;i++){
    rowMax[i]=-1;
  }
  IntegerVector rowMaxIdx(mat.nrow()*2-1);
  
  for(int i=0;i<mat.nrow();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.row(i)=mat.row(i);
  }
  for(int i=0;i<mat.nrow()-1;i++){
    for(int j=i+1;j<mat.nrow();j++){
      double shared=locusSUM(locusAND(mat.row(i),mat.row(j)),weight=weight);
      cacheShared[i][j]=shared;
      if(rowMax[i]<shared){
        rowMax[i]=shared;
        rowMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.nrow();
  while(nodeCount<mat.nrow()*2-1){
    if(verbose&&((nodeCount-mat.nrow())&1023)==1023){
      Rprintf("=");
    }
    
    double maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(rowMaxIdx(i))==-1){
          rowMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              double shared=cacheShared[cacheIdxi][cacheIdxj];
              if(rowMax[i]<shared){
                rowMax[i]=shared;
                rowMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<rowMax[i]){
          maxShared=rowMax[i];
          joinIdx1=i;
          joinIdx2=rowMaxIdx(i);
        }
      }
    }
    
    edge((nodeCount-mat.nrow())*2,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2,1)=joinIdx1;
    edge((nodeCount-mat.nrow())*2+1,0)=nodeCount;
    edge((nodeCount-mat.nrow())*2+1,1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerMatrix tmpMat(2,mat.ncol());
    tmpMat.row(0)=cacheMat.row(cacheIdx1);
    tmpMat.row(1)=cacheMat.row(cacheIdx2);
    IntegerVector stateAns=bitAncestral(tmpMat);
    
    for(int k=0;k<mat.ncol();k++){
      edge((nodeCount-mat.nrow())*2,3+k)=cacheMat(cacheIdx1,k);
      edge((nodeCount-mat.nrow())*2+1,3+k)=cacheMat(cacheIdx2,k);
    }
    
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.row(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        double shared=locusSUM(locusAND(cacheMat.row(tmp),cacheMat.row(cacheIdxAns)),weight=weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(rowMax[i]<shared){
          rowMax[i]=shared;
          rowMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

// [[Rcpp::export]]
IntegerVector radiusUnweightedParallel(int i,IntegerMatrix mat,int maxChar,int radius){
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerVector output(nbin+1);
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  int nNeighbor=0;
  for(int x=0;x<mat.ncol();x++){
    for(int y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
      }
    }
  }
  for(int j=0;j<mat.nrow();j++){
    int nshared=0;
    for(int k=0;k<nchar;k++){
      if(((mat(j,xidx[k])>>yidx[k])&1)==1){
        nshared++;
      }
    }
    if(nshared>=radius){
      int bin_x=j>>5;
      int bin_y=j&31;
      output(bin_x+1)|=(1<<bin_y);
      nNeighbor++;
    }
  }
  output(0)=nNeighbor;
  return output;
}

// [[Rcpp::export]]
IntegerVector radiusWeightedParallel(int i,IntegerMatrix mat,int maxChar,NumericVector weight,double radius){
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerVector output(nbin+1);
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  int nNeighbor=0;
  for(int x=0;x<mat.ncol();x++){
    for(int y=0;y<32;y++){
      if(((mat(i,x)>>y)&1)==1){
        xidx[nchar]=x;
        yidx[nchar]=y;
        nchar++;
      }
    }
  }
  for(int j=0;j<mat.nrow();j++){
    double nshared=0;
    for(int k=0;k<nchar;k++){
      if(((mat(j,xidx[k])>>yidx[k])&1)==1){
        nshared+=weight(xidx[k]*32+yidx[k]);
      }
    }
    if(nshared>=radius){
      int bin_x=j>>5;
      int bin_y=j&31;
      output(bin_x+1)|=(1<<bin_y);
      nNeighbor++;
    }
  }
  output(0)=nNeighbor;
  return output;
}
