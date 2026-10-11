#include <Rcpp.h>
#include <immintrin.h>
#include <chrono>
#include <random>
#include <string>

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

// // [[Rcpp::export]]
// std::vector<int> int2idx(IntegerVector ve) {
//   std::vector<int> output;
//   
//   const int n = ve.size();
//   for (int i = 0; i < n; ++i) {
//     unsigned int value = static_cast<unsigned int>(ve[i]);
//     int base = i * 32;
//     while (value != 0) {
//       unsigned int lowest = value & -value;
//       int j = 0;
//       unsigned int temp = lowest;
//       while (temp > 1) {
//         temp >>= 1;
//         j++;
//       }
//       output.push_back(base + j);
//       value &= value - 1;
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
std::vector<int> int2idx(IntegerVector ve) {
  std::vector<int> output;
  
  const int n = ve.size();
  for (int i = 0; i < n; ++i) {
    unsigned int value = static_cast<unsigned int>(ve[i]);
    int base = i * 32;
    while (value != 0) {
      int j = __builtin_ctz(value);
      output.push_back(base + j);
      value &= value - 1;
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
  float** cacheShared=(float**)(R_alloc(mat.nrow(),sizeof(float*)));
  for(int i=0;i<mat.nrow();i++){
    cacheShared[i]=(float*)(R_alloc(mat.nrow(),sizeof(float)));
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
    if(verbose&&((nodeCount-mat.nrow())&32767)==32767){
      Rprintf("\n");
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
  float** cacheShared=(float**)(R_alloc(mat.nrow(),sizeof(float*)));
  for(int i=0;i<mat.nrow();i++){
    cacheShared[i]=(float*)(R_alloc(mat.nrow(),sizeof(float)));
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
    if(verbose&&((nodeCount-mat.nrow())&32767)==32767){
      Rprintf("\n");
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

// // [[Rcpp::export]]
// IntegerVector edge2vector(IntegerMatrix edge,int maxIdx){
//   int i;
//   IntegerVector route(maxIdx);
//   route.fill(-1);
//   for(i=0;i<edge.nrow();i++){
//     route(edge(i,1)-1)=edge(i,0)-1;
//   }
//   return route;
// }

// [[Rcpp::export]]
IntegerVector edge2vector(IntegerMatrix edge,int len){
  //transform the edge into link vector to facilitate tip to root search
  IntegerVector linkVector(len);
  linkVector.fill(-1); // root remains the initial value, so it is linked to -1.
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    linkVector(desIdx)=ansIdx;
  }
  return linkVector;
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
    if(verbose&&(i&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&(i&32767)==32767){
      Rprintf("\n");
    }
    int idx[maxChar];
    int xidx[maxChar];
    int yidx[maxChar];
    int nchar=0;
    
    int cols = mat.ncol();
    for (int x = 0; x < cols; x++) {
      uint32_t val = mat(i, x);
      while (val) {
        int bit = __builtin_ctz(val);
        idx[nchar]=x*32+bit;
        xidx[nchar] = x;
        yidx[nchar] = bit;
        nchar++;
        val &= val - 1;
      }
    }
    
    // for(x=0;x<mat.ncol();x++){
    //   for(y=0;y<32;y++){
    //     if(((mat(i,x)>>y)&1)==1){
    //       idx[nchar]=x*32+y;
    //       xidx[nchar]=x;
    //       yidx[nchar]=y;
    //       nchar++;
    //     }
    //   }
    // }
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
    if(verbose&&(i&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&(i&32767)==32767){
      Rprintf("\n");
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
  int j,k;
  int nbin=((mat.nrow()-1)>>5)+1;
  IntegerVector mat_neighbor(nbin);
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  
  int cols = mat.ncol();
  for (int x = 0; x < cols; x++) {
    uint32_t val = mat(i, x);
    while (val) {
      int bit = __builtin_ctz(val);
      xidx[nchar] = x;
      yidx[nchar] = bit;
      nchar++;
      val &= val - 1;
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
    if(verbose&&(i&32767)==32767){
      Rprintf("\n");
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
  int neighbor_idx[n_cols];
  for(i=0;i<n_cols;i++){
    if(verbose&&(i&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&(i&32767)==32767){
      Rprintf("\n");
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
IntegerVector nodeDesBig(IntegerMatrix edge,int nnode,int ntip,std::string file_path,int cacheSize=8192,bool header=true){
  IntegerVector nodeSize(nnode);
  nodeSize.fill(0);
  int nbin=((ntip-1)>>5)+1;
  int** cacheMat=(int**)(R_alloc(cacheSize,sizeof(int*)));
  for(int i=0;i<cacheSize;i++){
    cacheMat[i]=(int*)(R_alloc(nbin,sizeof(int)));
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
  if(header){
    fwrite(dim,sizeof(int),2,file);    
  }

  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      nodeSize(ansIdx-ntip)++;
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      int ansCacheIdx=node2cache(ansIdx-ntip);
      if(ansCacheIdx==-1){
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return nodeSize;
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
      nodeSize(ansIdx-ntip)+=nodeSize(desIdx-ntip);
      if(node2cache(ansIdx-ntip)==-1){
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return nodeSize;
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
  IntegerVector output(nnode-1);
  for(int i=1;i<nnode;i++){
    output(nnode-1-i)=nodeSize(i);
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector findNeighborRangeUnweighted(int i,IntegerMatrix mat,int maxChar,IntegerMatrix range, int offset, int n=10,double p=0.5,int radius=0){
  IntegerVector output(range.ncol());
  int xidx[maxChar];
  int yidx[maxChar];
  int nchar=0;
  
  int cols = mat.ncol();
  for (int x = 0; x < cols; x++) {
    uint32_t val = mat(i, x);
    while (val) {
      int bit = __builtin_ctz(val);
      xidx[nchar] = x;
      yidx[nchar] = bit;
      nchar++;
      val &= val - 1;
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
  
  int cols = mat.ncol();
  for (int x = 0; x < cols; x++) {
    uint32_t val = mat(i, x);
    while (val) {
      int bit = __builtin_ctz(val);
      xidx[nchar] = x;
      yidx[nchar] = bit;
      nchar++;
      radiusSelf+=weight(x*32+bit);
      val &= val - 1;
    }
  }
  
  // for(int x=0;x<mat.ncol();x++){
  //   for(int y=0;y<32;y++){
  //     if(((mat(i,x)>>y)&1)==1){
  //       idx[nchar]=x*32+y;
  //       xidx[nchar]=x;
  //       yidx[nchar]=y;
  //       nchar++;
  //       radiusSelf+=weight(x*32+y);
  //     }
  //   }
  // }
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
    if(verbose&&((nodeCount-mat.nrow())&32767)==32767){
      Rprintf("\n");
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
    if(verbose&&((nodeCount-mat.nrow())&32767)==32767){
      Rprintf("\n");
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
IntegerVector radiusWeightedParallel(int i,IntegerMatrix mat,int maxChar,NumericVector weight,float radius){
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
    float nshared=0;
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

// [[Rcpp::export]]
NumericVector bsCollapse(IntegerMatrix edge,NumericVector edgeLen,NumericVector bs,double threshold=0.5){
  int nnode=bs.size()+1;
  int ntip=edgeLen.size()+1-nnode;
  NumericVector output=edgeLen;
  for(int i=0;i<edgeLen.size();i++){
    int idx=edge(i,1);
    if(idx>ntip&&bs(idx-1-ntip-1)<threshold){
      double len=output(i);
      output(i)=0;
      for(int j=i+1;j<edgeLen.size();j++){
        if(edge(j,0)==idx){
          output(j)+=len;
        }
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdge2(IntegerMatrix mat,NumericMatrix weight,IntegerVector tipSize,int minSize=0,float minRatio=2.0,bool verbose=false){
  // mat row:feature col:tip
  IntegerMatrix edge(3+mat.nrow(),mat.ncol()*2-1);
  IntegerVector node2cache(mat.ncol()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.ncol());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  float** cacheShared=(float**)(R_alloc(mat.ncol(),sizeof(float*)));
  for(int i=0;i<mat.ncol();i++){
    cacheShared[i]=(float*)(R_alloc(mat.ncol(),sizeof(float)));
  }
  float colMax[mat.ncol()*2-1];
  for(int i=0;i<mat.ncol()*2-1;i++){
    colMax[i]=-1;
  }
  IntegerVector colMaxIdx(mat.ncol()*2-1);
  IntegerVector nodeSize(mat.ncol()*2-1);
  for(int i=0;i<tipSize.size();i++){
    nodeSize(i)=tipSize(i);
  }
  
  for(int i=0;i<mat.ncol();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.column(i)=mat.column(i);
  }
  for(int i=0;i<mat.ncol()-1;i++){
    for(int j=i+1;j<mat.ncol();j++){
      float shared=countApoFloat(mat.column(i),mat.column(j),weight);
      cacheShared[i][j]=shared;
      if(colMax[i]<shared){
        colMax[i]=shared;
        colMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.ncol();
  while(nodeCount<mat.ncol()*2-1){
    if(verbose&&((nodeCount-mat.ncol())&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&((nodeCount-mat.ncol())&32767)==32767){
      Rprintf("\n");
    }
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(colMaxIdx(i))==-1){
          colMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(colMax[i]<shared){
                colMax[i]=shared;
                colMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<colMax[i]){
          maxShared=colMax[i];
          joinIdx1=i;
          joinIdx2=colMaxIdx(i);
        }
      }
    }
    
    edge(0,(nodeCount-mat.ncol())*2)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2)=joinIdx1;
    edge(0,(nodeCount-mat.ncol())*2+1)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2+1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerVector stateAns(mat.nrow());
    
    if(nodeSize(joinIdx1)<minSize&&nodeSize(joinIdx1)*minRatio<nodeSize(joinIdx2)){
      stateAns=cacheMat.column(cacheIdx2);
    }else if(nodeSize(joinIdx2)<minSize&&nodeSize(joinIdx2)*minRatio<nodeSize(joinIdx1)){
      stateAns=cacheMat.column(cacheIdx1);
    }else{
      for(int k=0;k<mat.nrow();k++){
        int tmp1=cacheMat(k,cacheIdx1);
        int tmp2=cacheMat(k,cacheIdx2);
        if(tmp1==-1){
          tmp1=tmp2;
        }else if(tmp2>=0&&tmp1!=tmp2){
          tmp1=0;
        }
        stateAns(k)=tmp1;
      }
    }
    nodeSize(nodeCount)=nodeSize(joinIdx1)+nodeSize(joinIdx2);
    
    for(int k=0;k<mat.nrow();k++){
      edge(3+k,(nodeCount-mat.ncol())*2)=cacheMat(k,cacheIdx1);
      edge(3+k,(nodeCount-mat.ncol())*2+1)=cacheMat(k,cacheIdx2);
      if(stateAns(k)==0&&cacheMat(k,cacheIdx1)>0){
        edge(2,(nodeCount-mat.ncol())*2)++;
      }
      if(stateAns(k)==0&&cacheMat(k,cacheIdx2)>0){
        edge(2,(nodeCount-mat.ncol())*2+1)++;
      }
    }
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.column(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=countApoFloat(cacheMat.column(tmp),cacheMat.column(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(colMax[i]<shared){
          colMax[i]=shared;
          colMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  
  for(int k=0;k<mat.nrow();k++){
    edge(3+k,edge.ncol()-1)=cacheMat(k,node2cache(nodeCount-1));
  }
  
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

// [[Rcpp::export]]
IntegerMatrix greedyJoinEdgeBit2(IntegerMatrix mat,NumericVector weight,IntegerVector tipSize,int minSize=0,float minRatio=2.0,bool verbose=false){
  // mat row:feature col:tip
  IntegerMatrix edge(3+mat.nrow(),mat.ncol()*2-1);
  IntegerVector node2cache(mat.ncol()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.ncol());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  float** cacheShared=(float**)(R_alloc(mat.ncol(),sizeof(float*)));
  for(int i=0;i<mat.ncol();i++){
    cacheShared[i]=(float*)(R_alloc(mat.ncol(),sizeof(float)));
  }
  float colMax[mat.ncol()*2-1];
  for(int i=0;i<mat.ncol()*2-1;i++){
    colMax[i]=-1;
  }
  IntegerVector colMaxIdx(mat.ncol()*2-1);
  IntegerVector nodeSize(mat.ncol()*2-1);
  for(int i=0;i<tipSize.size();i++){
    nodeSize(i)=tipSize(i);
  }
  
  for(int i=0;i<mat.ncol();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.column(i)=mat.column(i);
  }
  for(int i=0;i<mat.ncol()-1;i++){
    std::vector<int> idx=int2idx(mat.column(i));
    for(int j=i+1;j<mat.ncol();j++){
      float shared=0;
      for(unsigned int k=0;k<idx.size();k++){
        int xx=idx[k]>>5;
        int yy=idx[k]&31;
        if(((mat(xx,j)>>yy)&1)==1){
          shared+=weight(idx[k]);
        }
      }
      cacheShared[i][j]=shared;
      if(colMax[i]<shared){
        colMax[i]=shared;
        colMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.ncol();
  while(nodeCount<mat.ncol()*2-1){
    if(verbose&&((nodeCount-mat.ncol())&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&((nodeCount-mat.ncol())&32767)==32767){
      Rprintf("\n");
    }
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(colMaxIdx(i))==-1){
          colMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(colMax[i]<shared){
                colMax[i]=shared;
                colMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<colMax[i]){
          maxShared=colMax[i];
          joinIdx1=i;
          joinIdx2=colMaxIdx(i);
        }
      }
    }
    
    edge(0,(nodeCount-mat.ncol())*2)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2)=joinIdx1;
    edge(0,(nodeCount-mat.ncol())*2+1)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2+1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    
    IntegerVector stateAns;
    if(nodeSize(joinIdx1)<minSize&&nodeSize(joinIdx1)*minRatio<nodeSize(joinIdx2)){
      stateAns=cacheMat.column(cacheIdx2);
    }else if(nodeSize(joinIdx2)<minSize&&nodeSize(joinIdx2)*minRatio<nodeSize(joinIdx1)){
      stateAns=cacheMat.column(cacheIdx1);
    }else{
      stateAns=locusAND(cacheMat.column(cacheIdx1),cacheMat.column(cacheIdx2));
    }
    nodeSize(nodeCount)=nodeSize(joinIdx1)+nodeSize(joinIdx2);
    
    for(int k=0;k<mat.nrow();k++){
      edge(3+k,(nodeCount-mat.ncol())*2)=cacheMat(k,cacheIdx1);
      edge(3+k,(nodeCount-mat.ncol())*2+1)=cacheMat(k,cacheIdx2);
    }
    
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.column(cacheIdxAns)=stateAns;
    
    std::vector<int> idx=int2idx(cacheMat.column(cacheIdxAns));
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=0;
        for(unsigned int k=0;k<idx.size();k++){
          int xx=idx[k]>>5;
          int yy=idx[k]&31;
          if(((cacheMat(xx,tmp)>>yy)&1)==1){
            shared+=weight(idx[k]);
          }
        }
        // float shared=countApoBit(cacheMat.column(tmp),cacheMat.column(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(colMax[i]<shared){
          colMax[i]=shared;
          colMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  
  for(int k=0;k<mat.nrow();k++){
    edge(3+k,edge.ncol()-1)=cacheMat(k,node2cache(nodeCount-1));
  }
  if(verbose){
    Rprintf("\n");
  }
  return edge;
}

// // [[Rcpp::export]]
// IntegerVector nodeKeep2(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
//   // mat row:feature col:tip
//   int i,j,k=0;
//   int Nnode=node_mat.ncol();
//   IntegerVector state(Ntip+Nnode);
//   IntegerVector output(Ntip+Nnode);
//   state.fill(-1);
//   for(i=0;i<Ntip;i++){
//     state(i)=1;
//   }
//   for(i=edge.nrow()-1;i>=0;i--){
//     if(state(edge(i,1)-1)==0){
//       state(edge(i,0)-1)=0;
//     }else if(state(edge(i,1)-1)==1){
//       if(state(edge(i,0)-1)==-1){
//         state(edge(i,0)-1)=1;
//         for(j=0;j<outmat.ncol();j++){
//           if(mutationNested(outmat.column(j),node_mat.column(edge(i,0)-1-Ntip))<keep_cutoff){
//             state(edge(i,0)-1)=0;
//             break;
//           }
//         }
//       }
//       if(state(edge(i,0)-1)==0){
//         output(k)=edge(i,1);
//         k++;
//       }
//     }else if(state(edge(i,1)-1)==-1){
//       Rprintf("Error at row %d\n",i+1);
//       return state;
//     }
//   }
//   return output;
// }

// // [[Rcpp::export]]
// IntegerVector nodeKeepBit2(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
//   // mat row:feature col:tip
//   int i,j,k=0;
//   int Nnode=node_mat.ncol();
//   IntegerVector state(Ntip+Nnode);
//   IntegerVector output(Ntip+Nnode);
//   state.fill(-1);
//   for(i=0;i<Ntip;i++){
//     state(i)=1;
//   }
//   for(i=edge.nrow()-1;i>=0;i--){
//     if(state(edge(i,1)-1)==0){
//       state(edge(i,0)-1)=0;
//     }else if(state(edge(i,1)-1)==1){
//       if(state(edge(i,0)-1)==-1){
//         state(edge(i,0)-1)=1;
//         for(j=0;j<outmat.ncol();j++){
//           if(mutationNestedBit(outmat.column(j),node_mat.column(edge(i,0)-1-Ntip))<keep_cutoff){
//             state(edge(i,0)-1)=0;
//             break;
//           }
//         }
//       }
//       if(state(edge(i,0)-1)==0){
//         output(k)=edge(i,1);
//         k++;
//       }
//     }else if(state(edge(i,1)-1)==-1){
//       Rprintf("Error at row %d\n",i+1);
//       return state;
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
IntegerMatrix reOrderEdge2(IntegerMatrix edge,int root){
  // edge row:ans,des col:branch
  IntegerMatrix output(edge.nrow(),edge.ncol());
  IntegerVector ances_stack(edge.ncol()+1);
  ances_stack(0)=root;
  int stack_pointer=0;
  bool isfind;
  int i=0,j;
  while(stack_pointer>(-1)){
    isfind=false;
    for(j=0;j<edge.ncol();j++){
      if(edge(0,j)==ances_stack(stack_pointer)){
        output.column(i)=edge.column(j);
        i++;
        stack_pointer++;
        ances_stack(stack_pointer)=edge(1,j);
        edge(0,j)=-1;
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
std::vector<int> outmatFilterBit(IntegerMatrix mat,IntegerVector state,int cutoff=1){
  std::vector<int> xidx;
  std::vector<int> yidx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    for(int j=0;j<32;j++){
      if(((state(i)>>j)&1)==1){
        xidx.push_back(i);
        yidx.push_back(j);
      }
    }
  }
  for(int i=0;i<mat.ncol();i++){
    int nFound=0;
    for(unsigned int k=0;k<xidx.size();k++){
      if(((mat(xidx[k],i)>>yidx[k])&1)==0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(i);
    }
  }
  return output;
}

// [[Rcpp::export]]
int nodeFilterBit(IntegerMatrix mat,IntegerVector state,int cutoff=1){
  std::vector<int> xidx;
  std::vector<int> yidx;
  int output=1;
  for(int i=0;i<state.size();i++){
    for(int j=0;j<32;j++){
      if(((state(i)>>j)&1)==1){
        xidx.push_back(i);
        yidx.push_back(j);
      }
    }
  }
  for(int i=0;i<mat.ncol();i++){
    int nFound=0;
    for(unsigned int k=0;k<xidx.size();k++){
      if(((mat(xidx[k],i)>>yidx[k])&1)==0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output=0;
      break;
    }
  }
  return output;
}

// // [[Rcpp::export]]
// std::vector<int> nodeKeepBit2(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
//   // mat row:feature col:tip
//   int k=0;
//   int Nnode=node_mat.ncol();
//   IntegerVector state(Ntip+Nnode);
//   std::vector<int> output;
//   state.fill(-1);
//   for(int i=0;i<Ntip;i++){
//     state(i)=1;
//   }
//   for(int i=edge.nrow()-1;i>=0;i--){
//     if(state(edge(i,1)-1)==0){
//       state(edge(i,0)-1)=0;
//     }else if(state(edge(i,1)-1)==1){
//       if(state(edge(i,0)-1)==-1){
//         state(edge(i,0)-1)=nodeFilterBit(outmat,node_mat.column(edge(i,0)-1-Ntip),keep_cutoff);
//       }
//       if(state(edge(i,0)-1)==0){
//         output.push_back(edge(i,1));
//         k++;
//       }
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
std::vector<int> nodeKeepBit2(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
  // mat row:feature col:tip
  std::vector<int> output;
  int currentNode=Ntip;
  int retain=0;
  for(int i=0;i<edge.nrow();i++){
    if(edge(i,1)>Ntip){
      if(retain==1&&currentNode<=edge(i,0)){
        continue;
      }
      retain=nodeFilterBit(outmat,node_mat.column(edge(i,1)-1-Ntip),keep_cutoff);
      if(retain==1){
        currentNode=edge(i,1);
        output.push_back(currentNode);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> outmatFilterIndel(IntegerMatrix mat,IntegerVector state,int cutoff=1){
  std::vector<int> idx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    if(state(i)>0){
      idx.push_back(i);
    }
  }
  for(int i=0;i<mat.ncol();i++){
    int nFound=0;
    for(unsigned int k=0;k<idx.size();k++){
      if(mat(idx[k],i)!=state(idx[k])&&mat(idx[k],i)>=0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(i);
    }
  }
  return output;
}

// [[Rcpp::export]]
int nodeFilterIndel(IntegerMatrix mat,IntegerVector state,int cutoff=1){
  std::vector<int> idx;
  int output=1;
  for(int i=0;i<state.size();i++){
    if(state(i)>0){
      idx.push_back(i);
    }
  }
  for(int i=0;i<mat.ncol();i++){
    int nFound=0;
    for(unsigned int k=0;k<idx.size();k++){
      if(mat(idx[k],i)!=state(idx[k])&&mat(idx[k],i)>=0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output=0;
      break;
    }
  }
  return output;
}

// // [[Rcpp::export]]
// std::vector<int> nodeKeepIndel(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
//   // mat row:feature col:tip
//   int k=0;
//   int Nnode=node_mat.ncol();
//   IntegerVector state(Ntip+Nnode);
//   std::vector<int> output;
//   state.fill(-1);
//   for(int i=0;i<Ntip;i++){
//     state(i)=1;
//   }
//   for(int i=edge.nrow()-1;i>=0;i--){
//     if(state(edge(i,1)-1)==0){
//       state(edge(i,0)-1)=0;
//     }else if(state(edge(i,1)-1)==1){
//       if(state(edge(i,0)-1)==-1){
//         state(edge(i,0)-1)=nodeFilterIndel(outmat,node_mat.column(edge(i,0)-1-Ntip),keep_cutoff);
//       }
//       if(state(edge(i,0)-1)==0){
//         output.push_back(edge(i,1));
//         k++;
//       }
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
std::vector<int> nodeKeepIndel(IntegerMatrix edge,int Ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
  // mat row:feature col:tip
  std::vector<int> output;
  int currentNode=Ntip;
  int retain=0;
  for(int i=0;i<edge.nrow();i++){
    if(edge(i,1)>Ntip){
      if(retain==1&&currentNode<=edge(i,0)){
        continue;
      }
      retain=nodeFilterIndel(outmat,node_mat.column(edge(i,1)-1-Ntip),keep_cutoff);
      if(retain==1){
        currentNode=edge(i,1);
        output.push_back(currentNode);
      }
    }
  }
  return output;
}


// [[Rcpp::export]]
IntegerVector edgeSplit(IntegerMatrix edge,IntegerVector node){
  int k=0;
  bool isFree=true;
  IntegerVector category(edge.nrow());
  category.fill(-1);
  for(int i=0;i<edge.nrow();i++){
    if(isFree){
      if(edge(i,0)==node(k)){
        isFree=false;
        category(i)=k;
      }else{
      }
    }else{
      if(edge(i,0)>=node(k)){
        category(i)=k;
      }else{
        k++;
        if(k>=node.size()){
          break;
        }
        isFree=true;
      }
    }
  }
  return(category);
}

// [[Rcpp::export]]
IntegerVector tipMergeCross(SEXP pNeighborMat,int offset1,int offset2,int len,
                            IntegerVector nNeighbor,IntegerVector member,
                            IntegerVector transition,
                            double p=0.8,bool sqrt_nNeighbor=false,
                            int min_shared=0){
  XPtr<BigMatrix> xpMat(pNeighborMat);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int n_cols=xpMat->ncol();
  int n_rows=xpMat->nrow();
  IntegerVector output=Rcpp::clone(transition);
  
  IntegerVector neighborIdx(n_cols);
  IntegerVector xIdx(n_cols);
  IntegerVector yIdx(n_cols);
  for(int iTmp=0;iTmp<len&&iTmp+offset1<n_cols;iTmp++){
    int i=iTmp+offset1;
    int root1=0,root2=0;
    root1=transitionRoot(output,member(i));
    int nCount=0;
    for(int x=0;x<n_rows;x++){
      int mask=mat[i][x];
      for(int y=0;y<32&&nCount<nNeighbor(i);y++){
        if(((mask>>y)&1)==1){
          neighborIdx(nCount)=x*32+y;
          xIdx(nCount)=x;
          yIdx(nCount)=y;
          nCount++;
        }
      }
    }
    for(int j=0;j<nCount;j++){
      int currentIdx=neighborIdx(j);
      if(currentIdx>=offset2+len){
        break;
      }
      if(currentIdx<offset2||currentIdx<=i){
        continue;
      }
      root2=transitionRoot(output,member(currentIdx));
      if(root1==root2){
        continue;
      }
      int nshared=0;
      int mask=0;
      int ptr=-1;
      for(int k=0;k<nCount;k++){
        if(ptr!=xIdx(k)){
          mask=mat[currentIdx][xIdx(k)];
          ptr=xIdx(k);
        }
        if(((mask>>yIdx(k))&1)==1){
          nshared++;
        }
      }
      if(nshared<min_shared){
        continue;
      }
      double expected=0;
      if(sqrt_nNeighbor){
        expected=sqrt(nCount*nNeighbor(currentIdx))*p;
      }else{
        if(nCount<nNeighbor(currentIdx)){
          expected=nCount*p;
        }else{
          expected=nNeighbor(currentIdx)*p;
        }
      }
      if(nshared>expected){
        if(root1>root2){
          output(root1)=root2;
          root1=root2;
        }else{
          output(root2)=root1;
        }
      }
    }
  }
  for(int i=0;i<output.size();i++){
    output(i)=transitionRoot(output,output(i));
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector tipMergeSub(SEXP pNeighborMat,int offset,int len,
                          IntegerVector nNeighbor,
                          double p=0.8,bool sqrt_nNeighbor=false,
                          int min_shared=0){
  XPtr<BigMatrix> xpMat(pNeighborMat);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int n_cols=xpMat->ncol();
  int n_rows=xpMat->nrow();
  
  IntegerVector member(n_cols);
  member.fill(-1);
  IntegerVector transition(n_cols);
  transition.fill(-1);
  int nCategory=0;
  IntegerVector neighborIdx(n_cols);
  IntegerVector xIdx(n_cols);
  IntegerVector yIdx(n_cols);
  for(int iTmp=0;iTmp<len&&iTmp+offset<n_cols;iTmp++){
    int i=iTmp+offset;
    int root1=0,root2=0;
    if(member(i)==-1){
      root1=nCategory;
      member(i)=nCategory;
      transition(nCategory)=nCategory;
      nCategory++;
    }else{
      root1=transitionRoot(transition,member(i));
    }
    int nCount=0;
    for(int x=0;x<n_rows;x++){
      int mask=mat[i][x];
      for(int y=0;y<32&&nCount<nNeighbor(i);y++){
        if(((mask>>y)&1)==1){
          neighborIdx(nCount)=x*32+y;
          xIdx(nCount)=x;
          yIdx(nCount)=y;
          nCount++;
        }
      }
    }
    for(int j=0;j<nCount;j++){
      int currentIdx=neighborIdx(j);
      if(currentIdx>=offset+len){
        break;
      }
      if(i>=currentIdx||currentIdx<offset){
        continue;
      }
      root2=transitionRoot(transition,member(currentIdx));
      if(root1==root2){
        continue;
      }
      int nshared=0;
      int mask=0;
      int ptr=-1;
      for(int k=0;k<nCount;k++){
        if(ptr!=xIdx(k)){
          mask=mat[currentIdx][xIdx(k)];
          ptr=xIdx(k);
        }
        if(((mask>>yIdx(k))&1)==1){
          nshared++;
        }
      }
      if(nshared<min_shared){
        continue;
      }
      double expected=0;
      if(sqrt_nNeighbor){
        expected=sqrt(nCount*nNeighbor(currentIdx))*p;
      }else{
        if(nCount<nNeighbor(currentIdx)){
          expected=nCount*p;
        }else{
          expected=nNeighbor(currentIdx)*p;
        }
      }
      if(nshared>expected){
        if(root2==-1){
          member(currentIdx)=root1;
        }else{
          transition(root1)=root2;
          root1=root2;
        }
      }
    }
  }
  for(int i=0;i<member.size();i++){
    member(i)=transitionRoot(transition,member(i));
  }
  return member;
}

// [[Rcpp::export]]
IntegerVector combineTransition(IntegerVector transition1,IntegerVector transition2){
  IntegerVector output=Rcpp::clone(transition1);
  for(int i=0;i<output.size();i++){
    if(transition2(i)==i){
      continue;
    }
    int root1=transitionRoot(output,i);
    int root2=transitionRoot(output,transition2(i));
    if(root1>root2){
      output(root1)=root2;
    }else if(root1<root2){
      output(root2)=root1;
    }
  }
  for(int i=0;i<output.size();i++){
    output(i)=transitionRoot(output,output(i));
  }
  return output;
}

// // [[Rcpp::export]]
// std::vector<int> outmatFilterBit2(IntegerMatrix mat,IntegerVector rowIdx,IntegerVector state,int cutoff=1){
//   std::vector<int> xidx;
//   std::vector<int> yidx;
//   std::vector<int> output;
//   for(int i=0;i<state.size();i++){
//     for(int j=0;j<32;j++){
//       if(((state(i)>>j)&1)==1){
//         xidx.push_back(i);
//         yidx.push_back(j);
//       }
//     }
//   }
//   for(int i=0;i<rowIdx.size();i++){
//     int tmp=rowIdx(i);
//     int nFound=0;
//     for(unsigned int k=0;k<xidx.size();k++){
//       if(((mat(tmp,xidx[k])>>yidx[k])&1)==0){
//         nFound++;
//         if(nFound>=cutoff){
//           break;
//         }
//       }
//     }
//     if(nFound<cutoff){
//       output.push_back(tmp);
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
std::vector<int> outmatFilterBit2(IntegerMatrix mat,IntegerVector colIdx,IntegerVector state,int cutoff=1){
  std::vector<int> xidx;
  std::vector<int> yidx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    for(int j=0;j<32;j++){
      if(((state(i)>>j)&1)==1){
        xidx.push_back(i);
        yidx.push_back(j);
      }
    }
  }
  for(int i=0;i<colIdx.size();i++){
    int tmp=colIdx(i);
    int nFound=0;
    for(unsigned int k=0;k<xidx.size();k++){
      if(((mat(xidx[k],tmp)>>yidx[k])&1)==0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(tmp);
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> outmatFilterBitBig(SEXP bigmat,IntegerVector colIdx,IntegerVector state,int cutoff=1){
  XPtr<BigMatrix> xpMat(bigmat);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  std::vector<int> xidx;
  std::vector<int> yidx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    for(int j=0;j<32;j++){
      if(((state(i)>>j)&1)==1){
        xidx.push_back(i);
        yidx.push_back(j);
      }
    }
  }
  for(int i=0;i<colIdx.size();i++){
    int tmp=colIdx(i);
    int nFound=0;
    for(unsigned int k=0;k<xidx.size();k++){
      if(((mat[tmp][xidx[k]]>>yidx[k])&1)==0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(tmp);
    }
  }
  return output;
}

// // [[Rcpp::export]]
// std::vector<int> outmatFilterIndel2(IntegerMatrix mat,IntegerVector rowIdx,IntegerVector state,int cutoff=1){
//   std::vector<int> idx;
//   std::vector<int> output;
//   for(int i=0;i<state.size();i++){
//     if(state(i)>0){
//       idx.push_back(i);
//     }
//   }
//   for(int i=0;i<rowIdx.size();i++){
//     int nFound=0;
//     int tmp=rowIdx(i);
//     for(unsigned int k=0;k<idx.size();k++){
//       if(mat(tmp,idx[k])!=state(idx[k])&&mat(tmp,idx[k])>=0){
//         nFound++;
//         if(nFound>=cutoff){
//           break;
//         }
//       }
//     }
//     if(nFound<cutoff){
//       output.push_back(tmp);
//     }
//   }
//   return output;
// }

// [[Rcpp::export]]
std::vector<int> outmatFilterIndel2(IntegerMatrix mat,IntegerVector colIdx,IntegerVector state,int cutoff=1){
  std::vector<int> idx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    if(state(i)>0){
      idx.push_back(i);
    }
  }
  for(int i=0;i<colIdx.size();i++){
    int nFound=0;
    int tmp=colIdx(i);
    for(unsigned int k=0;k<idx.size();k++){
      if(mat(idx[k],tmp)!=state(idx[k])&&mat(idx[k],tmp)>=0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(tmp);
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> outmatFilterIndelBig(SEXP bigmat,IntegerVector colIdx,IntegerVector state,int cutoff=1){
  XPtr<BigMatrix> xpMat(bigmat);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  std::vector<int> idx;
  std::vector<int> output;
  for(int i=0;i<state.size();i++){
    if(state(i)>0){
      idx.push_back(i);
    }
  }
  for(int i=0;i<colIdx.size();i++){
    int nFound=0;
    int tmp=colIdx(i);
    for(unsigned int k=0;k<idx.size();k++){
      if(mat[tmp][idx[k]]!=state(idx[k])&&mat[tmp][idx[k]]>=0){
        nFound++;
        if(nFound>=cutoff){
          break;
        }
      }
    }
    if(nFound<cutoff){
      output.push_back(tmp);
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix subcols_bigmat(SEXP xp, Rcpp::IntegerVector cols) {
  Rcpp::XPtr<BigMatrix> xpMat(xp);
  BigMatrix* pMat = xpMat.get();
  
  const int nrow = pMat->nrow();
  const int ncol = pMat->ncol();
  const int ncol_new = cols.size();
  
  for (int j = 0; j < ncol_new; ++j) {
    if (cols[j] < 1 || cols[j] > ncol) {
      Rcpp::stop("col index out of range.");
    }
  }
  
  Rcpp::IntegerMatrix out(nrow, ncol_new);
  MatrixAccessor<int> acc(*pMat);
  
  for (int j = 0; j < ncol_new; ++j) {
    int* src = acc[cols[j] - 1];
    int* dst = out.begin() + static_cast<R_xlen_t>(j) * nrow;
    std::copy(src, src + nrow, dst);
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix selectedNodeDes(IntegerMatrix edge,int nnode,int ntip,IntegerVector node,int cacheSize=8192){
  int nbin=((ntip-1)>>5)+1;
  IntegerMatrix output(nbin,node.size());
  IntegerVector isOutput(nnode+ntip);
  isOutput.fill(-1);
  for(int i=0;i<node.size();i++){
    isOutput(node(i))=i;
    if(node(i)<ntip){
      output(node(i)>>5,i)=1<<(node(i)&31);
    }
  }
  // int** cacheMat=(int**)(malloc(cacheSize*sizeof(int*)));
  int** cacheMat=(int**)(R_alloc(cacheSize,sizeof(int*)));
  for(int i=0;i<cacheSize;i++){
    // cacheMat[i]=(int*)(malloc(nbin*sizeof(int)));
    cacheMat[i]=(int*)(R_alloc(nbin,sizeof(int)));
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
          return 0;
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
          return 0;
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
      int tmp=isOutput(desIdx);
      if(tmp!=-1){
        for(int x=0;x<nbin;x++){
          output(x,tmp)=cacheMat[desCacheIdx][x];
        }
      }
      cacheFree(desCacheIdx)=true;
      nUse--;
      for(int x=0;x<nbin;x++){
        cacheMat[desCacheIdx][x]=0;
      }
    }
  }
  int tmp=isOutput(ntip);
  if(tmp!=-1){
    for(int x=0;x<nbin;x++){
      output(x,tmp)=cacheMat[node2cache(0)][x];
    }
  }
  return output;
}

int countBitEvo(IntegerVector xx, IntegerVector yy) {
  int count=0;
  for(int i=0;i<yy.size();i++){
    int xorResult = xx(i) ^ yy(i);
    while (xorResult != 0) {
      xorResult = xorResult & (xorResult - 1);
      count++;
    }
  }
  return count;
}

// [[Rcpp::export]]
Rcpp::List bitAncestralEC(IntegerMatrix mat,IntegerMatrix edge,IntegerVector tipSize,int nnode,float minSize=10,float minRatio=2){
  int ntip=tipSize.size();
  IntegerVector nodeSize(nnode);
  nodeSize.fill(0);
  IntegerVector maxDes(nnode);
  maxDes.fill(0);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      nodeSize(ansIdx-ntip)+=tipSize(desIdx);
      if(maxDes(ansIdx-ntip)<tipSize(desIdx)){
        maxDes(ansIdx-ntip)=tipSize(desIdx);
      }
    }else{
      nodeSize(ansIdx-ntip)+=nodeSize(desIdx-ntip);
      if(maxDes(ansIdx-ntip)<nodeSize(desIdx-ntip)){
        maxDes(ansIdx-ntip)=nodeSize(desIdx-ntip);
      }
    }
  }
  IntegerMatrix nodeState(mat.nrow(),nnode);
  nodeState.fill(0xfffffff);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      if(minRatio*tipSize(desIdx)>=maxDes(ansIdx-ntip)||tipSize(desIdx)>=minSize){
        for(int j=0;j<mat.nrow();j++){
          nodeState(j,ansIdx-ntip)&=mat(j,desIdx);
        }
      }
    }else{
      if(minRatio*nodeSize(desIdx-ntip)>=maxDes(ansIdx-ntip)||nodeSize(desIdx-ntip)>=minSize){
        for(int j=0;j<mat.nrow();j++){
          nodeState(j,ansIdx-ntip)&=nodeState(j,desIdx-ntip);
        }
      }
    }
  }
  IntegerVector edgeLen(edge.nrow());
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      edgeLen(i)=countBitEvo(nodeState.column(ansIdx-ntip),mat.column(desIdx));
    }else{
      edgeLen(i)=countBitEvo(nodeState.column(ansIdx-ntip),nodeState.column(desIdx-ntip));
    }
  }
  return Rcpp::List::create(
    Named("nodeState")=nodeState,
    Named("edgeLen")=edgeLen
  );
}

int countCharEvo(IntegerVector xx, IntegerVector yy) {
  int count=0;
  for(int i=0;i<yy.size();i++){
    if(xx(i)>=0&&yy(i)>=0&&xx(i)!=yy(i)){
      count++;
    }
  }
  return count;
}

// [[Rcpp::export]]
Rcpp::List charAncestralEC(IntegerMatrix mat,IntegerMatrix edge,IntegerVector tipSize,int nnode,float minSize=0,float minRatio=0){
  int ntip=tipSize.size();
  IntegerVector nodeSize(nnode);
  nodeSize.fill(0);
  IntegerVector maxDes(nnode);
  maxDes.fill(0);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      nodeSize(ansIdx-ntip)+=tipSize(desIdx);
      if(maxDes(ansIdx-ntip)<tipSize(desIdx)){
        maxDes(ansIdx-ntip)=tipSize(desIdx);
      }
    }else{
      nodeSize(ansIdx-ntip)+=nodeSize(desIdx-ntip);
      if(maxDes(ansIdx-ntip)<nodeSize(desIdx-ntip)){
        maxDes(ansIdx-ntip)=nodeSize(desIdx-ntip);
      }
    }
  }
  IntegerMatrix nodeState(mat.nrow(),nnode);
  nodeState.fill(-1);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      if(minRatio*tipSize(desIdx)>=maxDes(ansIdx-ntip)||tipSize(desIdx)>=minSize){
        for(int j=0;j<mat.nrow();j++){
          if(nodeState(j,ansIdx-ntip)==-1){
            nodeState(j,ansIdx-ntip)=mat(j,desIdx);
          }else if(nodeState(j,ansIdx-ntip)!=mat(j,desIdx)&&mat(j,desIdx)!=-1){
            nodeState(j,ansIdx-ntip)=0;
          }
        }
      }
    }else{
      if(minRatio*nodeSize(desIdx-ntip)>=maxDes(ansIdx-ntip)||nodeSize(desIdx-ntip)>=minSize){
        for(int j=0;j<mat.nrow();j++){
          if(nodeState(j,ansIdx-ntip)==-1){
            nodeState(j,ansIdx-ntip)=nodeState(j,desIdx-ntip);
          }else if(nodeState(j,ansIdx-ntip)!=nodeState(j,desIdx-ntip)&&nodeState(j,desIdx-ntip)!=-1){
            nodeState(j,ansIdx-ntip)=0;
          }
        }
      }
    }
  }
  IntegerVector edgeLen(edge.nrow());
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      edgeLen(i)=countCharEvo(nodeState.column(ansIdx-ntip),mat.column(desIdx));
    }else{
      edgeLen(i)=countCharEvo(nodeState.column(ansIdx-ntip),nodeState.column(desIdx-ntip));
    }
  }
  return Rcpp::List::create(
    Named("nodeState")=nodeState,
    Named("edgeLen")=edgeLen
  );
}

// [[Rcpp::export]]
IntegerVector calPseudoSize(IntegerVector tipSize,IntegerVector tip2pseudo,int nPseudo) {
  IntegerVector output(nPseudo);
  output.fill(0);
  for(int i=0;i<tipSize.size();++i){
    int pseudoIdx=tip2pseudo[i]-1;
    output[pseudoIdx]+=tipSize[i];
  }
  return output;
}

float calcShared(std::vector<int> idx,IntegerVector seq,NumericVector weight){
  float shared=0;
  for(unsigned int k=0;k<idx.size();++k){
    int xx=idx[k]>>5;
    int yy=idx[k]&31;
    if((seq(xx)&(1<<yy))!=0){
      shared+=weight(idx[k]);
    }
  }
  return shared;
}

// // [[Rcpp::export]]
// IntegerMatrix calcAncestralBitReverse(
//     IntegerVector seq1,IntegerVector bistate1,
//     IntegerVector seq2,IntegerVector bistate2){
//   
//   int nbin=seq1.size();
//   IntegerVector seqRes(nbin);
//   IntegerVector bistateRes(nbin);
//   for(int i=0;i<nbin;++i){
//     const int s1 = seq1[i];
//     const int s2 = seq2[i];
//     const int bs1 = bistate1[i];
//     const int bs2 = bistate2[i];
//     
//     const int diff = s1 ^ s2; //s1 and s2 is different
//     const int conflict = diff & (bs1 | bs2); //different and one seq is bistate => res should be 0
//     
//     seqRes[i] = (s1 | s2) & ~conflict;
//     bistateRes[i] = (diff & ~(bs1 ^ bs2)) | (bs1 & bs2); //res is bistate if bs1=0, bs2=0,diff=1 or bs1=1, bs2=1
//   }
//   IntegerMatrix output(nbin,2);
//   output.column(0)=seqRes;
//   output.column(1)=bistateRes;
//   return(output);
// }

// [[Rcpp::export]]
IntegerMatrix calcAncestralBitReverse(
    IntegerVector seq1,IntegerVector bistate1,
    IntegerVector seq2,IntegerVector bistate2){
  
  int nbin=seq1.size();
  IntegerVector seqRes(nbin);
  IntegerVector bistateRes(nbin);
  
  size_t i = 0;
  for (; i + 8 <= nbin; i += 8) {
    __m256i s1_vec = _mm256_loadu_si256((__m256i*)&seq1[i]);
    __m256i s2_vec = _mm256_loadu_si256((__m256i*)&seq2[i]);
    __m256i bs1_vec = _mm256_loadu_si256((__m256i*)&bistate1[i]);
    __m256i bs2_vec = _mm256_loadu_si256((__m256i*)&bistate2[i]);
    
    __m256i diff_vec=_mm256_xor_si256(s1_vec,s2_vec);
    __m256i conflict_vec=_mm256_and_si256(_mm256_or_si256(bs1_vec,bs2_vec),diff_vec);
    __m256i seqRes_vec=_mm256_andnot_si256(conflict_vec,_mm256_or_si256(s1_vec,s2_vec));
    _mm256_storeu_si256((__m256i*)&seqRes[i], seqRes_vec);
    
    __m256i tmp1_vec=_mm256_andnot_si256(_mm256_xor_si256(bs1_vec,bs2_vec),diff_vec);
    __m256i tmp2_vec=_mm256_and_si256(bs1_vec,bs2_vec);
    __m256i bistateRes_vec=_mm256_or_si256(tmp1_vec,tmp2_vec);
    _mm256_storeu_si256((__m256i*)&bistateRes[i], bistateRes_vec);
  }
  
  for(;i<nbin;++i){
    const int s1 = seq1[i];
    const int s2 = seq2[i];
    const int bs1 = bistate1[i];
    const int bs2 = bistate2[i];
    
    const int diff = s1 ^ s2; //s1 and s2 is different
    const int conflict = diff & (bs1 | bs2); //different and one seq is bistate => res should be 0
    
    seqRes[i] = (s1 | s2) & ~conflict;
    bistateRes[i] = (diff & ~(bs1 ^ bs2)) | (bs1 & bs2); //res is bistate if bs1=0, bs2=0,diff=1 or bs1=1, bs2=1
  }  
  
  IntegerMatrix output(nbin,2);
  output.column(0)=seqRes;
  output.column(1)=bistateRes;
  return(output);
}


// // [[Rcpp::export]]
// IntegerMatrix calcAncestralBitReverse2(
//     IntegerVector seq1,IntegerVector bistate1,
//     IntegerVector seq2,IntegerVector bistate2){
//   
//   int nbin=seq1.size();
//   IntegerVector seqRes(nbin);
//   IntegerVector bistateRes(nbin);
//   for(int i=0;i<nbin;++i){
//     const int s1 = seq1[i];
//     const int s2 = seq2[i];
//     const int bs1 = bistate1[i];
//     const int bs2 = bistate2[i];
//     
//     seqRes[i]=(s1&s2)|(bs1&s2)|(bs2&s1);
//     bistateRes[i] = (s1^s2)&(~(bs1|bs2));
//   }
//   IntegerMatrix output(nbin,2);
//   output.column(0)=seqRes;
//   output.column(1)=bistateRes;
//   return(output);
// }

// [[Rcpp::export]]
IntegerMatrix calcAncestralBitReverse2(
    IntegerVector seq1,IntegerVector bistate1,
    IntegerVector seq2,IntegerVector bistate2){
  
  int nbin=seq1.size();
  IntegerVector seqRes(nbin);
  IntegerVector bistateRes(nbin);
  
  size_t i = 0;
  for (; i + 8 <= nbin; i += 8) {
    __m256i s1_vec = _mm256_loadu_si256((__m256i*)&seq1[i]);
    __m256i s2_vec = _mm256_loadu_si256((__m256i*)&seq2[i]);
    __m256i bs1_vec = _mm256_loadu_si256((__m256i*)&bistate1[i]);
    __m256i bs2_vec = _mm256_loadu_si256((__m256i*)&bistate2[i]);
    
    // (s1&s2)|(bs1&s2)|(bs2&s1)
    __m256i tmp1_vec=_mm256_and_si256(s1_vec,s2_vec);
    __m256i tmp2_vec=_mm256_and_si256(bs1_vec,s2_vec);
    __m256i tmp3_vec=_mm256_and_si256(bs2_vec,s1_vec);
    __m256i seqRes_vec=_mm256_or_si256(_mm256_or_si256(tmp1_vec,tmp2_vec),tmp3_vec);
    _mm256_storeu_si256((__m256i*)&seqRes[i], seqRes_vec);
    
    // (s1^s2)&(~(bs1|bs2))
    tmp1_vec=_mm256_xor_si256(s1_vec,s2_vec);
    tmp2_vec=_mm256_or_si256(bs1_vec,bs2_vec);
    __m256i bistateRes_vec=_mm256_andnot_si256(tmp2_vec,tmp1_vec);
    _mm256_storeu_si256((__m256i*)&bistateRes[i], bistateRes_vec);
  }
  
  for(;i<nbin;++i){
    const int s1 = seq1[i];
    const int s2 = seq2[i];
    const int bs1 = bistate1[i];
    const int bs2 = bistate2[i];
    
    seqRes[i]=(s1&s2)|(bs1&s2)|(bs2&s1);
    bistateRes[i] = (s1^s2)&(~(bs1|bs2));
  }
  
  IntegerMatrix output(nbin,2);
  output.column(0)=seqRes;
  output.column(1)=bistateRes;
  return(output);
}


// [[Rcpp::export]]
Rcpp::List greedyJoinEdgeBitReverse(
    IntegerMatrix mat,NumericVector weight,bool verbose=false, int ancestralType=0
){
  // mat row:feature col:tip
  int ntip=mat.ncol();
  int nbin=mat.nrow();
  IntegerMatrix edge(2,ntip*2-2);
  IntegerMatrix nodeState(nbin,ntip*2-1);
  IntegerMatrix nodeBistate(nbin,ntip*2-1);
  IntegerVector node2cache(ntip*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(ntip);
  IntegerMatrix cacheMat(nbin,ntip);
  IntegerMatrix cacheBistate(nbin,ntip);
  cacheBistate.fill(0);
  float** cacheShared=(float**)(R_alloc(ntip,sizeof(float*)));
  for(int i=0;i<ntip;i++){
    cacheShared[i]=(float*)(R_alloc(ntip,sizeof(float)));
  }
  // float colMax[mat.ncol()*2-1];
  float* colMax=(float*)(R_alloc(ntip*2-1,sizeof(float)));
  for(int i=0;i<ntip*2-1;i++){
    colMax[i]=-1;
  }
  IntegerVector colMaxIdx(ntip*2-1);
  
  for(int i=0;i<ntip;i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.column(i)=mat.column(i);
  }
  for(int i=0;i<ntip;i++){
    std::vector<int> idx=int2idx(mat.column(i));
    for(int j=i+1;j<mat.ncol();j++){
      float shared=calcShared(idx,mat.column(j),weight);
      cacheShared[i][j]=shared;
      if(colMax[i]<shared){
        colMax[i]=shared;
        colMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=ntip;
  while(nodeCount<ntip*2-1){
    if(verbose&&((nodeCount-ntip)&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&((nodeCount-ntip)&32767)==32767){
      Rprintf("\n");
    }
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(colMaxIdx(i))==-1){
          colMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(colMax[i]<shared){
                colMax[i]=shared;
                colMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<colMax[i]){
          maxShared=colMax[i];
          joinIdx1=i;
          joinIdx2=colMaxIdx(i);
        }
      }
    }
    
    edge(0,(nodeCount-mat.ncol())*2)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2)=joinIdx1;
    edge(0,(nodeCount-mat.ncol())*2+1)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2+1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerMatrix resAns;
    if(ancestralType==0){
      resAns=calcAncestralBitReverse(
        cacheMat.column(cacheIdx1),cacheBistate.column(cacheIdx1),
        cacheMat.column(cacheIdx2),cacheBistate.column(cacheIdx2)
      );
    }else{
      resAns=calcAncestralBitReverse2(
        cacheMat.column(cacheIdx1),cacheBistate.column(cacheIdx1),
        cacheMat.column(cacheIdx2),cacheBistate.column(cacheIdx2)
      );
    }
    
    IntegerVector stateAns=resAns.column(0);
    IntegerVector bistateAns=resAns.column(1);
    
    nodeState.column((nodeCount-ntip)*2)=cacheMat.column(cacheIdx1);
    nodeState.column((nodeCount-ntip)*2+1)=cacheMat.column(cacheIdx2);
    nodeBistate.column((nodeCount-ntip)*2)=cacheBistate.column(cacheIdx1);
    nodeBistate.column((nodeCount-ntip)*2+1)=cacheBistate.column(cacheIdx2);
    
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.column(cacheIdxAns)=stateAns;
    cacheBistate.column(cacheIdxAns)=bistateAns;
    
    std::vector<int> idx=int2idx(stateAns);
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=calcShared(idx,cacheMat.column(tmp),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(colMax[i]<shared){
          colMax[i]=shared;
          colMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  nodeState.column(ntip*2-2)=cacheMat.column(node2cache(nodeCount-1));
  nodeBistate.column(ntip*2-2)=cacheBistate.column(node2cache(nodeCount-1));
  
  if(verbose){
    Rprintf("\n");
  }
  
  for(int i=0;i<ntip*2-1;++i){
    for(int j=0;j<nbin;++j){
      nodeState(j,i)=nodeState(j,i)&(~nodeBistate(j,i));
    }
  }
  
  return Rcpp::List::create(
    Named("edge")=edge,
    Named("nodeState")=nodeState,
    Named("nodeBistate")=nodeBistate
  );
}

// [[Rcpp::export]]
IntegerVector calcDesStateBitReverse(
    IntegerVector seq1, IntegerVector seq2,IntegerVector bistate2) {
  
  IntegerVector output(seq1.size());
  output.fill(0);
  
  for (size_t i = 0; i < seq1.size(); ++i) {
    int s1 = seq1(i);
    int s2 = seq2(i);
    int bs2=bistate2(i);
    output(i)=(s2&~bs2)|(s1&bs2);
  }
  return output;
}

// [[Rcpp::export]]
Rcpp::List ancestralBitReverse(IntegerMatrix mat,IntegerMatrix edge,int ntip,int nnode){
  IntegerMatrix nodeState(mat.nrow(),nnode);
  nodeState.fill(0xffffffff);
  IntegerMatrix nodeBistate(mat.nrow(),nnode);
  nodeBistate.fill(0xffffffff);
  IntegerMatrix nodeConflict(mat.nrow(),nnode);
  nodeConflict.fill(0);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      for(int j=0;j<mat.nrow();j++){
        int s1=nodeState(j,ansIdx-ntip);
        int s2=mat(j,desIdx);
        int bs1=nodeBistate(j,ansIdx-ntip);
        int bs2=0x0;
        nodeState(j,ansIdx-ntip)=s1&s2;
        nodeBistate(j,ansIdx-ntip)&=bs2;
        nodeConflict(j,ansIdx-ntip)|=(s1^s2)&~(bs1|bs2);
      }
    }else{
      for(int j=0;j<mat.nrow();j++){
        nodeState(j,desIdx-ntip)|=nodeConflict(j,desIdx-ntip);
        nodeBistate(j,desIdx-ntip)|=nodeConflict(j,desIdx-ntip);
        
        int s1=nodeState(j,ansIdx-ntip);
        int s2=nodeState(j,desIdx-ntip);
        int bs1=nodeBistate(j,ansIdx-ntip);
        int bs2=nodeBistate(j,desIdx-ntip);
        nodeState(j,ansIdx-ntip)=s1&s2;
        nodeBistate(j,ansIdx-ntip)&=bs2;
        nodeConflict(j,ansIdx-ntip)|=(s1^s2)&~(bs1|bs2);
      }
    }
  }
  //update root
  for(int j=0;j<mat.nrow();j++){
    nodeState(j,0)&=~(nodeConflict(j,0)|nodeBistate(j,0));
  }
  
  IntegerVector edgeLen(edge.nrow());
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      int count=0;
      for(int j=0;j<mat.nrow();++j){
        edgeLen(i)=countBitEvo(nodeState.column(ansIdx-ntip),mat.column(desIdx));
      }
    }else{
      IntegerVector tmp=calcDesStateBitReverse(
        nodeState.column(ansIdx-ntip),
        nodeState.column(desIdx-ntip),nodeBistate.column(desIdx-ntip)
      );
      nodeState.column(desIdx-ntip)=tmp;
      edgeLen(i)=countBitEvo(nodeState.column(ansIdx-ntip),tmp);
    }
  }
  return Rcpp::List::create(
    Named("nodeState")=nodeState,
    Named("edgeLen")=edgeLen
  );
}

// [[Rcpp::export]]
IntegerVector encodeATGC8(const std::string& s){
  size_t n = s.size();
  if (n == 0) return IntegerVector(0);
  
  size_t out_len = (n + 3)/4;
  IntegerVector out(out_len);
  
  auto encode = [](char c) -> unsigned char {
    switch (c) {
    case 'A': case 'a': return 0x01;
    case 'T': case 't': return 0x02;
    case 'C': case 'c': return 0x04;
    case 'G': case 'g': return 0x08;
    case 'W': case 'w': return 0x03; //AT
    case 'M': case 'm': return 0x05; //AC
    case 'R': case 'r': return 0x09; //AG
    case 'Y': case 'y': return 0x06; //TC
    case 'K': case 'k': return 0x0A; //TG
    case 'S': case 's': return 0x0C; //CG
    case 'H': case 'h': return 0x07; //ATC
    case 'D': case 'd': return 0x0B; //ATG
    case 'V': case 'v': return 0x0D; //ACG
    case 'B': case 'b': return 0x0E; //TCG
    case 'N': case 'n': return 0x0F; //ATCG
    default: return 0x10;
    }
  };
  
  for (size_t i = 0; i < out_len; ++i) {
    int v = 0;
    for (int j = 0; j < 4; ++j) {
      size_t pos = i * 4 + j;
      if (pos < n) {
        unsigned char code = encode(s[pos]);
        v |= (code << (j * 8));
      }
    }
    out[i] = v;
  }
  return out;
}

// [[Rcpp::export]]
std::string decodeATGC8(const IntegerVector& encoded, bool randomize_ambiguous = false) {
  size_t n = encoded.size();
  if (n == 0) return "";
  
  size_t original_len = n * 4;
  std::string result;
  result.reserve(original_len);
  
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  
  auto resolve_ambiguous = [&](char ambig_base) -> char {
    if (!randomize_ambiguous) return ambig_base;
    
    std::uniform_int_distribution<int> dist_two(0, 1);
    std::uniform_int_distribution<int> dist_three(0, 2);
    std::uniform_int_distribution<int> dist_four(0, 3);
    
    switch (ambig_base) {
    case 'W': 
      return dist_two(generator) == 0 ? 'A' : 'T';
    case 'M':
      return dist_two(generator) == 0 ? 'A' : 'C';
    case 'R':
      return dist_two(generator) == 0 ? 'A' : 'G';
    case 'Y':
      return dist_two(generator) == 0 ? 'C' : 'T';
    case 'K':
      return dist_two(generator) == 0 ? 'G' : 'T';
    case 'S':
      return dist_two(generator) == 0 ? 'C' : 'G';
    case 'H': {
        int r = dist_three(generator);
        if (r == 0) return 'A';
        else if (r == 1) return 'C';
        else return 'T';
      }
    case 'D': {
      int r = dist_three(generator);
      if (r == 0) return 'A';
      else if (r == 1) return 'G';
      else return 'T';
    }
    case 'V': {
      int r = dist_three(generator);
      if (r == 0) return 'A';
      else if (r == 1) return 'C';
      else return 'G';
    }
    case 'B': {
      int r = dist_three(generator);
      if (r == 0) return 'C';
      else if (r == 1) return 'G';
      else return 'T';
    }
    case 'N': {
      int r = dist_four(generator);
      switch (r) {
      case 0: return 'A';
      case 1: return 'C';
      case 2: return 'G';
      default: return 'T';
      }
    }
    default:
      return ambig_base;
    }
  };
  
  auto decode = [&](unsigned int code) -> char {
    char base;
    switch (code) {
    case 0x01: base = 'A'; break;
    case 0x02: base = 'T'; break;
    case 0x04: base = 'C'; break;
    case 0x08: base = 'G'; break;
    case 0x03: base = 'W'; break;
    case 0x05: base = 'M'; break;
    case 0x09: base = 'R'; break;
    case 0x06: base = 'Y'; break;
    case 0x0A: base = 'K'; break;
    case 0x0C: base = 'S'; break;
    case 0x07: base = 'H'; break;
    case 0x0B: base = 'D'; break;
    case 0x0D: base = 'V'; break;
    case 0x0E: base = 'B'; break;
    case 0x0F: base = 'N'; break;
    default: base = '-'; break;
    }
    return resolve_ambiguous(base);
  };
  
  for (size_t i = 0; i < n; ++i) {
    unsigned int value = encoded[i];
    
    for (int j = 0; j < 4; ++j) {
      unsigned int code = (value >> (j * 8)) & 0xFF;
      result += decode(code);
    }
  }
  return result;
}

// [[Rcpp::export]]
IntegerVector countATCG(String s) {
  std::string str = s;  
  IntegerVector cnt = IntegerVector(4);
  cnt.fill(0);
  for (char ch : str) {
    switch (std::tolower(ch)) {
    case 'a': cnt(0)++; break;
    case 't': cnt(1)++; break;
    case 'c': cnt(2)++; break;
    case 'g': cnt(3)++; break;
    default:  break;
    }
  }
  return cnt;
}

// [[Rcpp::export]]
IntegerVector calcAncestorState(IntegerVector seq1, IntegerVector seq2) {
  
  IntegerVector ancestor(seq1.size());
  ancestor.fill(0);
  
  for (size_t i = 0; i < seq1.size(); ++i) {
    unsigned int value1 = seq1[i];
    unsigned int value2 = seq2[i];
    unsigned int ancestor_value = 0;
    
    for (size_t j = 0; j < 4; ++j) {
      unsigned char nuc1 = (value1 >> (j * 8)) & 0xFF;
      unsigned char nuc2 = (value2 >> (j * 8)) & 0xFF;
      unsigned char ancesTmp = 0;
      
      unsigned char intersection = nuc1 & nuc2;
      if (intersection != 0) {
        ancesTmp = intersection;
      } else {
        ancesTmp = nuc1 | nuc2;
      }
      
      ancestor_value |= (ancesTmp << (j * 8));
    }
    
    ancestor[i] = ancestor_value;
  }
  return ancestor;
}


// [[Rcpp::export]]
float seqSimilarity(IntegerVector seq1, IntegerVector seq2, NumericVector weight) {
  float total_similarity = 0.0;
  for (size_t i = 0; i < seq1.size(); ++i) {
    unsigned int value=seq1[i]&seq2[i];
    for (size_t j = 0; j < 4; ++j) {
      unsigned char nuc = (value >> (j * 8)) & 0x0F;
      size_t weight_index = i * 4 + j;
      if (nuc != 0x00) {
        total_similarity += (float)(weight[weight_index]);
      }
    }
  }
  return total_similarity;
}

// [[Rcpp::export]]
Rcpp::List greedyJoinEdgeDNA(IntegerMatrix mat,NumericVector weight,bool verbose=false){
  // mat row:feature col:tip
  IntegerMatrix edge(2,mat.ncol()*2-2);
  IntegerMatrix nodeState(mat.nrow(),mat.ncol()*2-1);
  IntegerVector node2cache(mat.ncol()*2-1);
  node2cache.fill(-1);
  IntegerVector cache2node(mat.ncol());
  IntegerMatrix cacheMat(mat.nrow(),mat.ncol());
  float** cacheShared=(float**)(R_alloc(mat.ncol(),sizeof(float*)));
  for(int i=0;i<mat.ncol();i++){
    cacheShared[i]=(float*)(R_alloc(mat.ncol(),sizeof(float)));
  }
  float colMax[mat.ncol()*2-1];
  for(int i=0;i<mat.ncol()*2-1;i++){
    colMax[i]=-1;
  }
  IntegerVector colMaxIdx(mat.ncol()*2-1);
  
  for(int i=0;i<mat.ncol();i++){
    node2cache(i)=i;
    cache2node(i)=i;
    cacheMat.column(i)=mat.column(i);
  }
  for(int i=0;i<mat.ncol()-1;i++){
    for(int j=i+1;j<mat.ncol();j++){
      float shared=seqSimilarity(mat.column(i),mat.column(j),weight);
      cacheShared[i][j]=shared;
      if(colMax[i]<shared){
        colMax[i]=shared;
        colMaxIdx(i)=j;
      }
    }
  }
  
  int nodeCount=mat.ncol();
  while(nodeCount<mat.ncol()*2-1){
    if(verbose&&((nodeCount-mat.ncol())&1023)==1023){
      Rprintf("=");
    }
    if(verbose&&((nodeCount-mat.ncol())&32767)==32767){
      Rprintf("\n");
    }
    float maxShared=-1;
    int joinIdx1=0,joinIdx2=0;
    for(int i=0;i<nodeCount-1;i++){
      int cacheIdxi=node2cache(i);
      if(cacheIdxi!=-1){
        if(node2cache(colMaxIdx(i))==-1){
          colMax[i]=-1;
          for(int j=i+1;j<nodeCount;j++){
            int cacheIdxj=node2cache(j);
            if(cacheIdxj!=-1){
              float shared=cacheShared[cacheIdxi][cacheIdxj];
              if(colMax[i]<shared){
                colMax[i]=shared;
                colMaxIdx(i)=j;
              }
            }
          }
        }
        if(maxShared<colMax[i]){
          maxShared=colMax[i];
          joinIdx1=i;
          joinIdx2=colMaxIdx(i);
        }
      }
    }
    
    edge(0,(nodeCount-mat.ncol())*2)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2)=joinIdx1;
    edge(0,(nodeCount-mat.ncol())*2+1)=nodeCount;
    edge(1,(nodeCount-mat.ncol())*2+1)=joinIdx2;
    
    int cacheIdx1=node2cache(joinIdx1);
    int cacheIdx2=node2cache(joinIdx2);
    IntegerVector stateAns=calcAncestorState(cacheMat.column(cacheIdx1),cacheMat.column(cacheIdx2));
    
    nodeState.column((nodeCount-mat.ncol())*2)=cacheMat.column(cacheIdx1);
    nodeState.column((nodeCount-mat.ncol())*2+1)=cacheMat.column(cacheIdx2);
    
    int cacheIdxAns=cacheIdx1;
    node2cache(nodeCount)=cacheIdxAns;
    cache2node(cacheIdxAns)=nodeCount;
    node2cache(joinIdx1)=-1;
    node2cache(joinIdx2)=-1;
    cacheMat.column(cacheIdxAns)=stateAns;
    
    for(int i=0;i<nodeCount;i++){
      int tmp=node2cache(i);
      if(tmp!=-1){
        float shared=seqSimilarity(cacheMat.column(tmp),cacheMat.column(cacheIdxAns),weight);
        cacheShared[tmp][cacheIdxAns]=shared;
        if(colMax[i]<shared){
          colMax[i]=shared;
          colMaxIdx(i)=nodeCount;
        }
      }
    }
    nodeCount++;
  }
  nodeState.column(nodeState.ncol()-1)=cacheMat.column(node2cache(nodeCount-1));
  
  if(verbose){
    Rprintf("\n");
  }
  
  return Rcpp::List::create(
    Named("edge")=edge,
    Named("nodeState")=nodeState
  );
}

// [[Rcpp::export]]
IntegerVector calcDesState(IntegerVector seq1, IntegerVector seq2) {
  IntegerVector output(seq1.size());
  output.fill(0);
  
  for (size_t i = 0; i < seq1.size(); ++i) {
    unsigned int value1 = seq1[i];
    unsigned int value2 = seq2[i];
    unsigned int output_value = 0;
    
    for (size_t j = 0; j < 4; ++j) {
      unsigned char nuc1 = (value1 >> (j * 8)) & 0xFF;
      unsigned char nuc2 = (value2 >> (j * 8)) & 0xFF;
      unsigned char outputTmp = 0;
      
      unsigned char intersection = nuc1 & nuc2;
      if (intersection == 0x0) {
        outputTmp = nuc2;
      } else {
        outputTmp = intersection;
      }
      
      output_value |= (outputTmp << (j * 8));
    }
    
    output[i] = output_value;
  }
  return output;
}

std::vector<int> countNucleoEvo(IntegerVector seq1, IntegerVector seq2) {
  std::vector<int> output;
  for (size_t i = 0; i < seq1.size(); ++i) {
    unsigned int value1 = seq1[i];
    unsigned int value2 = seq2[i];
    
    for (size_t j = 0; j < 4; ++j) {
      unsigned char nuc1 = (value1 >> (j * 8)) & 0xFF;
      unsigned char nuc2 = (value2 >> (j * 8)) & 0xFF;
      
      if (nuc1 != 0 || nuc2 != 0) {
        if ((nuc1 & nuc2) == 0) {
          output.push_back(i*4+j);
        }
      }
    }
  }
  return output;
}

IntegerVector getSeqState(IntegerVector seqIntersect, IntegerVector seqUnion) {
  IntegerVector output(seqIntersect.size());
  output.fill(0);
  
  for (size_t i = 0; i < seqIntersect.size(); i++) {
    unsigned int nuc_intersect = seqIntersect(i);
    unsigned int nuc_union = seqUnion(i);
    unsigned int output_value = 0;
    
    for (size_t j = 0; j < 4; j++) {
      unsigned char value_intersect = (nuc_intersect >> (j * 8)) & 0xFF;
      unsigned char value_union = (nuc_union >> (j * 8)) & 0xFF;
      unsigned char res = 0x00;
      
      if (value_intersect == 0) {
        res = value_union;
      } else {
        res = value_intersect;
      }
      
      output_value |= (res << (j * 8));
    }
    
    output(i) = output_value;
  }
  return output;
}


// [[Rcpp::export]]
Rcpp::List nucleoAncestral(IntegerMatrix mat,IntegerMatrix edge,int ntip,int nnode){
  IntegerMatrix nodeState(mat.nrow(),nnode);
  IntegerMatrix nodeStateIntersect(mat.nrow(),nnode);
  nodeStateIntersect.fill(0xffffffff);
  IntegerMatrix nodeStateUnion(mat.nrow(),nnode);
  nodeStateUnion.fill(0);
  
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      for(int j=0;j<mat.nrow();j++){
        nodeStateIntersect(j,ansIdx-ntip)&=mat(j,desIdx);
        nodeStateUnion(j,ansIdx-ntip)|=mat(j,desIdx);
      }
    }else{
      nodeState.column(desIdx-ntip)=getSeqState(
        nodeStateIntersect.column(desIdx-ntip),nodeStateUnion.column(desIdx-ntip)
      );
      for(int j=0;j<mat.nrow();j++){
        nodeStateIntersect(j,ansIdx-ntip)&=nodeState(j,desIdx-ntip);
        nodeStateUnion(j,ansIdx-ntip)|=nodeState(j,desIdx-ntip);
      }
    }
  }
  //update root
  nodeState.column(0)=getSeqState(
    nodeStateIntersect.column(0),nodeStateUnion.column(0)
  );
  IntegerMatrix nodeState_intermediate=clone(nodeState);
  
  IntegerVector edgeLen(edge.nrow());
  IntegerVector nEvo(mat.nrow()*4);
  nEvo.fill(0);
  IntegerVector nEvo2(mat.nrow()*4);
  nEvo2.fill(0);
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      std::vector<int> siteIdx=countNucleoEvo(nodeState.column(ansIdx-ntip),mat.column(desIdx));
      edgeLen(i)=siteIdx.size();
      for(size_t j=0;j<siteIdx.size();++j){
        nEvo[siteIdx[j]]++;
      }
    }else{
      IntegerVector tmp=calcDesState(nodeState.column(ansIdx-ntip),nodeState.column(desIdx-ntip));
      nodeState.column(desIdx-ntip)=tmp;
      std::vector<int> siteIdx=countNucleoEvo(nodeState.column(ansIdx-ntip),tmp);
      edgeLen(i)=siteIdx.size();
      for(size_t j=0;j<siteIdx.size();++j){
        nEvo[siteIdx[j]]++;
        int xx=siteIdx[j]/4;
        int yy=siteIdx[j]&3;
        unsigned char nuc1=(nodeState(xx,ansIdx-ntip)>>(yy*8))&0xFF;
        unsigned char nuc2=(nodeState(xx,desIdx-ntip)>>(yy*8))&0xFF;
        if(nuc1!=0x10&&nuc2!=0x10){
          nEvo2[siteIdx[j]]++;
        }
      }
    }
  }
  return Rcpp::List::create(
    Named("nodeState")=nodeState,
    Named("nodeState_intermediate")=nodeState_intermediate,
    Named("edgeLen")=edgeLen,
    Named("nEvo")=nEvo,
    Named("nEvo2")=nEvo2
  );
}

// [[Rcpp::export]]
IntegerVector radiusWeightedDNAParallel(int i,IntegerMatrix mat,NumericVector weight,float radius){
  int nbin=(mat.ncol()+31)>>5;
  IntegerVector output(nbin+1);
  int nNeighbor=0;
  for(int j=0;j<mat.ncol();j++){
    float shared=seqSimilarity(mat.column(i),mat.column(j),weight);
    if(shared>=radius){
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
std::vector<int> outmatFilterDNA(IntegerMatrix mat,IntegerVector colIdx,IntegerVector state,int cutoff=1){
  std::vector<int> output;
  int monoCount=0;
  int polyCount=0;
  IntegerVector monoIdx(state.size()*4);
  IntegerVector polyIdx(state.size()*4);
  for(int j=0;j<state.size();++j){
    int value=state[j];
    for(int k=0;k<4;++k){
      unsigned char nuc=(value>>(k*8))&0xFF;
      if(nuc==0x01||nuc==0x02||nuc==0x04||nuc==0x08){
        monoIdx[monoCount]=j*4+k;
        monoCount++;
      }else if(nuc>0&&nuc!=0x10){
        polyIdx[polyCount]=j*4+k;
        polyCount++;
      }
    }
  }
  for(int i=0;i<colIdx.size();++i){
    int tmp=colIdx(i);
    int nFound=0;
    for(int j=0;j<monoCount;++j){
      int idx=monoIdx[j];
      int xx=idx/4;
      int yy=idx&3;
      unsigned char nuc1=(state[xx]>>(yy*8))&0x0F;
      unsigned char nuc2=(mat(xx,tmp)>>(yy*8))&0x0F;
      if((nuc1&nuc2)==0){
        nFound++;
      }
      if(nFound-polyCount>=cutoff){
        break;
      }
    }
    if(nFound-polyCount>=cutoff){
      continue;
    }
    for(int j=0;j<polyCount;++j){
      int idx=polyIdx[j];
      int xx=idx/4;
      int yy=idx&3;
      unsigned char nuc1=(state[xx]>>(yy*8))&0x0F;
      unsigned char nuc2=(mat(xx,tmp)>>(yy*8))&0x0F;
      if((nuc1&nuc2)!=0){
        nFound--;
      }
    }
    if(nFound<cutoff){
      output.push_back(tmp);
    }
  }
  return output;
}

// [[Rcpp::export]]
int nodeFilterDNA(IntegerMatrix mat,IntegerVector state,int cutoff=1){
  int monoCount=0;
  int polyCount=0;
  IntegerVector monoIdx(state.size()*4);
  IntegerVector polyIdx(state.size()*4);
  for(int j=0;j<state.size();++j){
    int value=state[j];
    for(int k=0;k<4;++k){
      unsigned char nuc=(value>>(k*8))&0xFF;
      if(nuc==0x01||nuc==0x02||nuc==0x04||nuc==0x08){
        monoIdx[monoCount]=j*4+k;
        monoCount++;
      }else if(nuc>0&&nuc!=0x10){
        polyIdx[polyCount]=j*4+k;
        polyCount++;
      }
    }
  }
  for(int i=0;i<mat.ncol();++i){
    int nFound=0;
    for(int j=0;j<monoCount;++j){
      int idx=monoIdx[j];
      int xx=idx/4;
      int yy=idx&3;
      unsigned char nuc1=(state[xx]>>(yy*8))&0x0F;
      unsigned char nuc2=(mat(xx,i)>>(yy*8))&0x0F;
      if((nuc1&nuc2)==0){
        nFound++;
      }
      if(nFound-polyCount>=cutoff){
        break;
      }
    }
    if(nFound-polyCount>=cutoff){
      continue;
    }
    for(int j=0;j<polyCount;++j){
      int idx=polyIdx[j];
      int xx=idx/4;
      int yy=idx&3;
      unsigned char nuc1=(state[xx]>>(yy*8))&0x0F;
      unsigned char nuc2=(mat(xx,i)>>(yy*8))&0x0F;
      if((nuc1&nuc2)!=0){
        nFound--;
      }
    }
    if(nFound<cutoff){
      return 0;
    }
  }
  return 1;
}

// [[Rcpp::export]]
std::vector<int> nodeKeepDNA(IntegerMatrix edge,int ntip,IntegerMatrix node_mat,IntegerMatrix outmat,int keep_cutoff){
  // mat row:feature col:tip
  std::vector<int> output;
  int currentNode=ntip;
  int retain=0;
  for(int i=0;i<edge.nrow();i++){
    if(edge(i,1)>ntip){
      if(retain==1&&currentNode<=edge(i,0)){
        continue;
      }
      retain=nodeFilterDNA(outmat,node_mat.column(edge(i,1)-1-ntip),keep_cutoff);
      if(retain==1){
        currentNode=edge(i,1);
        output.push_back(currentNode);
      }
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix findRoot(IntegerMatrix edge,IntegerMatrix nodeState,IntegerMatrix mat,IntegerVector outIdx){
  
  int nNodes=nodeState.ncol();
  int nPositions=nodeState.nrow();
  int nSamples=outIdx.size();
  
  IntegerMatrix edgeState(nPositions,nNodes);
  edgeState.fill(0);
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    for(int j=0;j<nPositions;++j)
      edgeState(j,desIdx)=nodeState(j,ansIdx)|nodeState(j,desIdx);
  }
  IntegerVector rootFound(nSamples);
  IntegerVector nShared(nSamples);
  for(int sIdx=0;sIdx<nSamples;++sIdx){
    int maxShared=-1;
    int maxNode=-1;
    for(int i=0;i<edge.nrow();++i){
      int desIdx=edge(i,1)-1;
      int count=0;
      for(int j=0;j<nPositions;++j){
        unsigned int value=mat(j,outIdx[sIdx])&edgeState(j,desIdx);
        for(int k=0;k<4;++k){
          unsigned char nuc=(value>>(k*8))&0x0F;
          if(nuc!=0x00){
            count++;
          }
        }
      }
      if(count>maxShared){
        maxShared=count;
        maxNode=desIdx;
      }
    }
    rootFound[sIdx]=maxNode+1;
    nShared[sIdx]=maxShared;
  }
  IntegerMatrix output(nSamples,2);
  output.column(0)=rootFound;
  output.column(1)=nShared;
  return output;
}

// [[Rcpp::export]]
IntegerMatrix findRootAVX2(IntegerMatrix edge, IntegerMatrix nodeState, IntegerMatrix mat,IntegerVector outIdx) {
  
  int nNodes = nodeState.ncol();
  int nPositions = nodeState.nrow();
  int nSamples = outIdx.size();
  
  IntegerMatrix edgeState(nPositions, nNodes);
  edgeState.fill(0);
  
  for(int i = 0; i < edge.nrow(); ++i) {
    int ansIdx = edge(i, 0) - 1;
    int desIdx = edge(i, 1) - 1;
    for(int j = 0; j < nPositions; ++j) {
      edgeState(j, desIdx) = nodeState(j, ansIdx) | nodeState(j, desIdx);
    }
    // const int* ansData = &nodeState(0, ansIdx);
    // const int* desData = &nodeState(0, desIdx);
    // int j = 0;
    // for(; j + 7 < nPositions; j += 8) {
    //   __m256i ansVec = _mm256_loadu_si256((__m256i*)(ansData + j));
    //   __m256i desVec = _mm256_loadu_si256((__m256i*)(desData + j));
    //   __m256i orResult = _mm256_or_si256(ansVec, desVec);
    //   _mm256_storeu_si256((__m256i*)(&edgeState(j, desIdx)), orResult);
    // }
    // for(; j < nPositions; ++j) {
    //   edgeState(j, desIdx) = nodeState(j, ansIdx) | nodeState(j, desIdx);
    // }
  }
  
  IntegerVector rootFound(nSamples);
  IntegerVector nShared(nSamples);
  
  auto computeShared = [&](int outIdx, int desIdx) -> int {
    int count = 0;
    
    const int* outData = &mat(0, outIdx);
    const int* edgeData = &edgeState(0, desIdx);
    
    int i = 0;
    for(; i + 7 < nPositions; i += 8) {
      __m256i outVec = _mm256_loadu_si256((__m256i*)(outData + i));
      __m256i edgeVec = _mm256_loadu_si256((__m256i*)(edgeData + i));
      __m256i andResult = _mm256_and_si256(outVec, edgeVec);
      
      unsigned int results[8];
      _mm256_storeu_si256((__m256i*)results, andResult);
      
      for(int k = 0; k < 8; ++k) {
        if(results[k] != 0) {
          unsigned int value = results[k];
          for(int byte = 0; byte < 4; ++byte) {
            unsigned char nuc = (value >> (byte * 8)) & 0x0F;
            if(nuc != 0x00) {
              count++;
            }
          }
        }
      }
    }
    
    for(; i < nPositions; ++i) {
      unsigned int value = outData[i] & edgeData[i];
      if(value != 0) {
        for(int byte = 0; byte < 4; ++byte) {
          unsigned char nuc = (value >> (byte * 8)) & 0x0F;
          if(nuc != 0x00) {
            count++;
          }
        }
      }
    }
    
    return count;
  };
  
  for(int sIdx = 0; sIdx < nSamples; ++sIdx) {
    int maxShared = -1;
    int maxNode = -1;
    
    for(int i = 0; i < edge.nrow(); ++i) {
      int desIdx = edge(i, 1) - 1;
      
      int count = computeShared(outIdx[sIdx], desIdx);
      
      if(count > maxShared) {
        maxShared = count;
        maxNode = desIdx;
      }
    }
    
    rootFound[sIdx] = maxNode + 1;
    nShared[sIdx] = maxShared;
  }
  
  IntegerMatrix output(nSamples, 2);
  output.column(0) = rootFound;
  output.column(1) = nShared;
  return output;
}

// [[Rcpp::export]]
IntegerMatrix findRootAVX2Big(IntegerMatrix edge, IntegerMatrix nodeState, 
                              SEXP matSEXP, IntegerVector outIdx) {
  XPtr<BigMatrix> xpMat(matSEXP);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  
  int nNodes = nodeState.ncol();
  int nPositions = nodeState.nrow();
  int nSamples = outIdx.size();
  
  IntegerMatrix edgeState(nPositions, nNodes);
  edgeState.fill(0);
  
  for(int i = 0; i < edge.nrow(); ++i) {
    int ansIdx = edge(i, 0) - 1;
    int desIdx = edge(i, 1) - 1;
    for(int j = 0; j < nPositions; ++j) {
      edgeState(j, desIdx) = nodeState(j, ansIdx) | nodeState(j, desIdx);
    }
  }
  
  IntegerVector rootFound(nSamples);
  IntegerVector nShared(nSamples);
  std::vector<int> tempBuffer(nPositions);
  
  auto computeShared = [&](int outIdx, int desIdx) -> int {
    int count = 0;
    
    const int* edgeData = &edgeState(0, desIdx);
    const int* outData = mat[outIdx];
    int i = 0;
    for(; i + 7 < nPositions; i += 8) {
      __m256i outVec = _mm256_loadu_si256((__m256i*)(outData + i));
      __m256i edgeVec = _mm256_loadu_si256((__m256i*)(edgeData + i));
      __m256i andResult = _mm256_and_si256(outVec, edgeVec);
      
      unsigned int results[8];
      _mm256_storeu_si256((__m256i*)results, andResult);
      
      for(int k = 0; k < 8; ++k) {
        if(results[k] != 0) {
          unsigned int value = results[k];
          for(int byte = 0; byte < 4; ++byte) {
            unsigned char nuc = (value >> (byte * 8)) & 0x0F;
            if(nuc != 0x00) {
              count++;
            }
          }
        }
      }
    }
    
    for(; i < nPositions; ++i) {
      unsigned int value = mat[outIdx][i] & edgeData[i];
      if(value != 0) {
        for(int byte = 0; byte < 4; ++byte) {
          unsigned char nuc = (value >> (byte * 8)) & 0x0F;
          if(nuc != 0x00) {
            count++;
          }
        }
      }
    }
    
    return count;
  };
  
  for(int sIdx = 0; sIdx < nSamples; ++sIdx) {
    int maxShared = -1;
    int maxNode = -1;
    
    for(int i = 0; i < edge.nrow(); ++i) {
      int desIdx = edge(i, 1) - 1;
      
      int count = computeShared(outIdx[sIdx], desIdx);
      
      if(count > maxShared) {
        maxShared = count;
        maxNode = desIdx;
      }
    }
    
    rootFound[sIdx] = maxNode + 1;
    nShared[sIdx] = maxShared;
  }
  
  IntegerMatrix output(nSamples, 2);
  output.column(0) = rootFound;
  output.column(1) = nShared;
  return output;
}

// [[Rcpp::export]]
std::vector<int> breakUnrootedPhy(IntegerMatrix edge,int ntip,int nnode,IntegerVector breakNode){
  IntegerVector linkVector=edge2vector(edge,ntip+nnode);
  IntegerVector nodeState(ntip+nnode);
  nodeState.fill(1);
  nodeState[ntip]=0;
  for(int i=0;i<breakNode.size();++i){
    int nodePass=breakNode[i]-1;
    while(nodePass!=ntip){
      nodePass=linkVector[nodePass];
      nodeState[nodePass]=0;
    }
  }
  std::vector<int> output;
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(nodeState[ansIdx]==0&&nodeState[desIdx]==1){
      output.push_back(desIdx+1);
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix getEdgeState(IntegerMatrix edge,IntegerMatrix nodeState,IntegerMatrix outmat){
  
  int nNodes=nodeState.ncol();
  int nPositions=nodeState.nrow();
  
  IntegerMatrix edgeState(nPositions,nNodes);
  edgeState.fill(0);
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    for(int j=0;j<nPositions;++j)
      edgeState(j,desIdx)=nodeState(j,ansIdx)|nodeState(j,desIdx);
  }
  return edgeState;
}

// [[Rcpp::export]]
IntegerVector locusOR(IntegerVector ve1,IntegerVector ve2){
  IntegerVector output(ve1.size());
  for(int i=0;i<ve1.size();++i){
    output(i)=ve1(i)|ve2(i);
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector zipMat(IntegerVector ve){
  IntegerVector output((ve.size()-1)/2+1);
  for(int i=0;i<output.size();++i){
    int valueLower=ve(i*2);
    int valueUpper;
    if(i*2+1<ve.size()){
      valueUpper=ve(i*2+1);
    }else{
      valueUpper=0;
    }
    int valueZipped=0;
    for(int j=0;j<4;++j){
      int nuc=(valueLower>>(j*8))&0x0f;
      valueZipped|=(nuc<<(j*4));
    }
    for(int j=0;j<4;++j){
      int nuc=(valueUpper>>(j*8))&0x0f;
      valueZipped|=(nuc<<((j+4)*4));
    }
    output(i)=valueZipped;
  }
  return output;
}

// [[Rcpp::export]]
IntegerVector idx2weight(IntegerVector idx,int len){
  IntegerVector output(len);
  output.fill(0);
  for(int i=0;i<idx.size();++i){
    ++output(idx(i)-1);
  }
  return output;
}


// [[Rcpp::export]]
Rcpp::IntegerVector calcDesStateBitReverse_AVX2(
    Rcpp::IntegerVector seq1, 
    Rcpp::IntegerVector seq2,
    Rcpp::IntegerVector bistate2) {
  
  size_t n = seq1.size();
  Rcpp::IntegerVector output(n);
  if (n == 0) return output;
  
  size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    __m256i s1_vec = _mm256_loadu_si256((__m256i*)&seq1[i]);
    __m256i s2_vec = _mm256_loadu_si256((__m256i*)&seq2[i]);
    __m256i bs2_vec = _mm256_loadu_si256((__m256i*)&bistate2[i]);
    
    // output = (s2 & ~bs2) | (s1 & bs2)
    // output = (bs2 & s1) | (~bs2 & s2)
    // output = s2 ^ (bs2 & (s1 ^ s2))
    
    __m256i xor_vec = _mm256_xor_si256(s1_vec, s2_vec); // s1 ^ s2
    __m256i masked = _mm256_and_si256(xor_vec, bs2_vec); // bs2 & (s1 ^ s2)
    __m256i result = _mm256_xor_si256(s2_vec, masked);   // s2 ^ (bs2 & (s1 ^ s2))
    _mm256_storeu_si256((__m256i*)&output[i], result);
  }
  for (; i < n; ++i) {
    int s1 = seq1[i];
    int s2 = seq2[i];
    int bs2 = bistate2[i];
    output[i] = (s2 & ~bs2) | (s1 & bs2);
  }
  return output;
}

static inline void processTipNode_AVX2(
    int* nodeState_ans,
    const int* mat_col,
    int* nodeBistate_ans,
    int* nodeConflict_ans,
    size_t nrow) {
  
  size_t i = 0;
  __m256i zero = _mm256_setzero_si256();
  __m256i all_ones = _mm256_set1_epi32(-1);
  
  for (; i + 8 <= nrow; i += 8) {
    __m256i s1 = _mm256_loadu_si256((__m256i*)&nodeState_ans[i]);
    __m256i s2 = _mm256_loadu_si256((__m256i*)&mat_col[i]);
    __m256i bs1 = _mm256_loadu_si256((__m256i*)&nodeBistate_ans[i]);
    __m256i conf = _mm256_loadu_si256((__m256i*)&nodeConflict_ans[i]);
    
    // 1. nodeState = s1 & s2
    __m256i new_state = _mm256_and_si256(s1, s2);
    
    // 2. nodeBistate &= 0 
    __m256i new_bistate = zero;
    
    // 3. nodeConflict |= (s1 ^ s2) & ~bs1
    __m256i xor_result = _mm256_xor_si256(s1, s2);
    __m256i not_bs1 = _mm256_andnot_si256(bs1, all_ones);
    __m256i conflict = _mm256_and_si256(xor_result, not_bs1);
    __m256i new_conf = _mm256_or_si256(conf, conflict);
    
    _mm256_storeu_si256((__m256i*)&nodeState_ans[i], new_state);
    _mm256_storeu_si256((__m256i*)&nodeBistate_ans[i], new_bistate);
    _mm256_storeu_si256((__m256i*)&nodeConflict_ans[i], new_conf);
  }
  
  for (; i < nrow; ++i) {
    int s1 = nodeState_ans[i];
    int s2 = mat_col[i];
    int bs1 = nodeBistate_ans[i];
    
    nodeState_ans[i] = s1 & s2;
    nodeBistate_ans[i] = 0;  // & 0
    nodeConflict_ans[i] |= (s1 ^ s2) & ~bs1;
  }
}

static inline void processInternalNode_AVX2(
    int* nodeState_ans, int* nodeState_des,
    int* nodeBistate_ans, int* nodeBistate_des,
    int* nodeConflict_ans, int* nodeConflict_des,
    size_t nrow) {
  
  size_t i = 0;
  const __m256i all_ones = _mm256_set1_epi32(-1);
  
  for (; i + 8 <= nrow; i += 8) {
    __m256i s1 = _mm256_loadu_si256((const __m256i*)&nodeState_ans[i]);
    __m256i s2 = _mm256_loadu_si256((const __m256i*)&nodeState_des[i]);
    __m256i bs1 = _mm256_loadu_si256((const __m256i*)&nodeBistate_ans[i]);
    __m256i bs2 = _mm256_loadu_si256((const __m256i*)&nodeBistate_des[i]);
    __m256i conf_ans = _mm256_loadu_si256((const __m256i*)&nodeConflict_ans[i]);
    __m256i conf_des = _mm256_loadu_si256((const __m256i*)&nodeConflict_des[i]);
    
    
    __m256i new_s2 = _mm256_or_si256(s2, conf_des);
    __m256i new_bs2 = _mm256_or_si256(bs2, conf_des);
    __m256i new_s1 = _mm256_and_si256(s1, new_s2);
    __m256i new_bs1 = _mm256_and_si256(bs1, new_bs2);
    __m256i xor_result = _mm256_xor_si256(s1, new_s2);
    __m256i bs_or = _mm256_or_si256(bs1, new_bs2);
    __m256i not_bs_or = _mm256_andnot_si256(bs_or, all_ones);
    __m256i conflict = _mm256_and_si256(xor_result, not_bs_or);
    __m256i new_conf_ans = _mm256_or_si256(conf_ans, conflict);
    
    _mm256_storeu_si256((__m256i*)&nodeState_ans[i], new_s1);
    _mm256_storeu_si256((__m256i*)&nodeBistate_ans[i], new_bs1);
    _mm256_storeu_si256((__m256i*)&nodeState_des[i], new_s2);
    _mm256_storeu_si256((__m256i*)&nodeBistate_des[i], new_bs2);
    _mm256_storeu_si256((__m256i*)&nodeConflict_ans[i], new_conf_ans);
  }
  
  for (; i < nrow; ++i) {
    nodeState_des[i] |= nodeConflict_des[i];
    nodeBistate_des[i] |= nodeConflict_des[i];
    
    int s1 = nodeState_ans[i];
    int s2 = nodeState_des[i];
    int bs1 = nodeBistate_ans[i];
    int bs2 = nodeBistate_des[i];
    
    nodeState_ans[i] = s1 & s2;
    nodeBistate_ans[i] &= bs2;
    nodeConflict_ans[i] |= (s1 ^ s2) & ~(bs1 | bs2);
  }
}

static inline __m256i popcount_swar_avx2(__m256i v) {
  const __m256i mask1 = _mm256_set1_epi32(0x55555555);
  const __m256i mask2 = _mm256_set1_epi32(0x33333333);
  const __m256i mask4 = _mm256_set1_epi32(0x0f0f0f0f);
  const __m256i mask8 = _mm256_set1_epi32(0x00ff00ff);
  const __m256i mask16 = _mm256_set1_epi32(0x0000ffff);
  
  __m256i shifted = _mm256_srli_epi32(v, 1);
  shifted = _mm256_and_si256(shifted, mask1);
  v = _mm256_sub_epi32(v, shifted);
  
  __m256i temp = _mm256_and_si256(v, mask2);
  shifted = _mm256_srli_epi32(v, 2);
  shifted = _mm256_and_si256(shifted, mask2);
  v = _mm256_add_epi32(temp, shifted);
  
  shifted = _mm256_srli_epi32(v, 4);
  v = _mm256_add_epi32(v, shifted);
  v = _mm256_and_si256(v, mask4);
  
  shifted = _mm256_srli_epi32(v, 8);
  v = _mm256_add_epi32(v, shifted);
  v = _mm256_and_si256(v, mask8);
  
  shifted = _mm256_srli_epi32(v, 16);
  v = _mm256_add_epi32(v, shifted);
  v = _mm256_and_si256(v, mask16);
  
  return v;
}

static inline int popcount_scalar(int x) {
  int count = 0;
  while (x != 0) {
    x = x & (x - 1);
    count++;
  }
  return count;
}

static inline int count_hamming_avx2(const int* a, const int* b, size_t n) {
  size_t i = 0;
  int count = 0;
  __m256i sum_vec = _mm256_setzero_si256();
  
  for (; i + 8 <= n; i += 8) {
    __m256i va = _mm256_loadu_si256((const __m256i*)&a[i]);
    __m256i vb = _mm256_loadu_si256((const __m256i*)&b[i]);
    __m256i vxor = _mm256_xor_si256(va, vb);
    __m256i pop = popcount_swar_avx2(vxor);
    sum_vec = _mm256_add_epi32(sum_vec, pop);
  }
  
  int sum_array[8];
  _mm256_storeu_si256((__m256i*)sum_array, sum_vec);
  count = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3] +
    sum_array[4] + sum_array[5] + sum_array[6] + sum_array[7];
  
  for (; i < n; ++i) {
    count += popcount_scalar(a[i] ^ b[i]);
  }
  
  return count;
}

// [[Rcpp::export]]
Rcpp::List ancestralBitReverse_AVX2(IntegerMatrix mat,IntegerMatrix edge,int ntip,int nnode){
  int nrow = mat.nrow();
  IntegerMatrix nodeState(mat.nrow(),nnode);
  nodeState.fill(0xffffffff);
  IntegerMatrix nodeBistate(mat.nrow(),nnode);
  nodeBistate.fill(0xffffffff);
  IntegerMatrix nodeConflict(mat.nrow(),nnode);
  nodeConflict.fill(0);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      int* nodeState_ans = &nodeState(0, ansIdx - ntip);
      const int* mat_col = &mat(0, desIdx);
      int* nodeBistate_ans = &nodeBistate(0, ansIdx - ntip);
      int* nodeConflict_ans = &nodeConflict(0, ansIdx - ntip);
      
      processTipNode_AVX2(
        nodeState_ans, mat_col,
        nodeBistate_ans, nodeConflict_ans,
        nrow
      );
    }else{
      int* nodeState_ans = &nodeState(0, ansIdx - ntip);
      int* nodeState_des = &nodeState(0, desIdx - ntip);
      int* nodeBistate_ans = &nodeBistate(0, ansIdx - ntip);
      int* nodeBistate_des = &nodeBistate(0, desIdx - ntip);
      int* nodeConflict_ans = &nodeConflict(0, ansIdx - ntip);
      int* nodeConflict_des = &nodeConflict(0, desIdx - ntip);
      
      processInternalNode_AVX2(
        nodeState_ans, nodeState_des,
        nodeBistate_ans, nodeBistate_des,
        nodeConflict_ans, nodeConflict_des,
        nrow
      );
    }
  }
  //update root
  for(int j=0;j<mat.nrow();j++){
    nodeState(j,0)&=~(nodeConflict(j,0)|nodeBistate(j,0));
  }
  
  IntegerVector edgeLen(edge.nrow());
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      const int* state_col = &nodeState(0, ansIdx - ntip);
      const int* mat_col = &mat(0, desIdx);
      edgeLen(i) = count_hamming_avx2(state_col, mat_col, nrow);
    }else{
      IntegerVector tmp = calcDesStateBitReverse_AVX2(
        nodeState.column(ansIdx - ntip),
        nodeState.column(desIdx - ntip),
        nodeBistate.column(desIdx - ntip)
      );
      nodeState.column(desIdx-ntip)=tmp;
      const int* state_col = &nodeState(0, ansIdx - ntip);
      const int* tmp_ptr = tmp.begin();
      edgeLen(i) = count_hamming_avx2(state_col, tmp_ptr, nrow);
    }
  }
  return Rcpp::List::create(
    Named("nodeState")=nodeState,
    Named("edgeLen")=edgeLen
  );
}


int bitCompare(int* xx,IntegerVector yy){
  //if identical return 0, if xx>yy return 1, if xx<yy return 2
  for(int i=0;i<yy.size();i++){
    int tmp1=xx[i]&(~yy(i));
    int tmp2=(~xx[i])&yy(i);
    if(tmp1!=0){
      return(1);
    }else if(tmp2!=0){
      return(2);
    } 
  }
  return(0);
}

// [[Rcpp::export]]
Rcpp::List findMaxNested(IntegerVector clade,IntegerMatrix edge,int ntip,int nnode){
  //for a given set of tips (clade),find largest nodes in the tree (edge) within it
  IntegerVector nodeState(ntip+nnode);
  //first assume all nodes within the clade
  nodeState.fill(1);
  IntegerVector nodeSize(ntip+nnode);
  nodeSize.fill(0);
  for(int i=edge.nrow()-1;i>=0;i--){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(desIdx<ntip){
      //if desIdx is a tip
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      if(((clade(x)>>y)&1)==0){
        //if desIdx is not in clade, mark both ans and des as outside the clade(nodeState=0)
        nodeState(ansIdx)=0;
        nodeState(desIdx)=0;
      }else{
        nodeSize(desIdx)=1;
        nodeSize(ansIdx)+=1;
      }
    }else{
      //if desIdx is a node
      if(nodeState(desIdx)==0){
        // if desIdx is not in clade, mark ans as not in clade(nodeState=0)
        nodeState(ansIdx)=0;
      }else{
        nodeSize(ansIdx)+=nodeSize(desIdx);
      }
    } 
  }
  //after loop, nodeState of all nodes outside clade marked as 0 
  
  //search for edge that ans outside clade, des in clade, des is a largest node in clade
  std::vector<int> output;
  std::vector<int> outputSize;
  for(int i=0;i<edge.nrow();i++){
    int ansIdx=edge(i,0)-1;
    int desIdx=edge(i,1)-1;
    if(nodeState(ansIdx)==0&&nodeState(desIdx)==1){
      output.push_back(desIdx+1);
      outputSize.push_back(nodeSize(desIdx));
    }
  }
  return Rcpp::List::create(Named("node")=output,Named("size")=outputSize);
}

IntegerVector countDifferentBits(int* xx, IntegerVector yy) {
  // count different bits (count(0)) and unique bits in xx (count(1))
  IntegerVector count(2);
  count.fill(0);
  for(int i=0;i<yy.size();i++){
    int xorResult = xx[i] ^ yy(i);
    while (xorResult != 0) {
      xorResult = xorResult & (xorResult - 1);
      count(0)++;
    }
    xorResult = xx[i] & (~yy(i));
    while (xorResult != 0) {
      xorResult = xorResult & (xorResult - 1);
      count(1)++;
    }
  }
  return count;
}

// [[Rcpp::export]]
Rcpp::List approximateCompare(SEXP bitClade,int ntip,
                              IntegerMatrix edge1,
                              IntegerMatrix edge2,int nnode2,
                              int cacheSize=8192,float threshold=0.1){
  XPtr<BigMatrix> xpMat(bitClade);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int nnode1=xpMat->ncol()+1;
  int nbin=xpMat->nrow();
  IntegerMatrix cacheMat(nbin,cacheSize);
  cacheMat.fill(0);
  LogicalVector cacheFree(cacheSize);
  cacheFree.fill(true);
  int nUse=0;
  IntegerVector node2cache(nnode2);
  node2cache.fill(-1);
  int cacheIdx=0;
  
  IntegerVector outputDis1(nnode1-1);
  outputDis1.fill(ntip*2);
  IntegerVector outputNode1(nnode1-1);
  outputNode1.fill(0);
  IntegerVector outputDis2(nnode2);
  outputDis2.fill(ntip*2);
  IntegerVector outputNode2(nnode2);
  outputNode2.fill(0);
  
  //create search idx for edge1
  IntegerVector linkVector=edge2vector(edge1,ntip+nnode1);
  // IntegerVector linkVector(ntip+nnode1);
  // for(int i=0;i<edge1.nrow();i++){
  //   int ansIdx=edge1(i,0)-1;
  //   int desIdx=edge1(i,1)-1;
  //   linkVector(desIdx)=ansIdx;
  // }
  //the node of edge1 where node of edge2 start to search
  IntegerMatrix startNode(nnode2,2);
  startNode.fill(-1);
  //the sizes of nodes in edge2
  IntegerVector nodeSize2(ntip+nnode2);
  nodeSize2.fill(0);
  
  for(int i=edge2.nrow()-1;i>=0;i--){
    int ansIdx=edge2(i,0)-1;
    int desIdx=edge2(i,1)-1;
    if(desIdx<ntip){ // des is a tip
      nodeSize2(desIdx)=1;
      nodeSize2(ansIdx)++;
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      int ansCacheIdx=node2cache(ansIdx-ntip); // find ans in cache
      if(ansCacheIdx==-1){ //if ans is not in cache,store it in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        ansCacheIdx=cacheIdx;
        node2cache(ansIdx-ntip)=cacheIdx; // link node to cache
        cacheFree(cacheIdx)=false; //mark cache as used
        nUse++;
      }
      cacheMat(x,ansCacheIdx)|=(1<<y); // update the member of ans
      // when ansIdx start searching, it can begin from the ancestor of desIdx in edge1
      if(startNode(ansIdx-ntip,0)==-1){
        startNode(ansIdx-ntip,0)=linkVector(desIdx);
      }else{
        startNode(ansIdx-ntip,1)=linkVector(desIdx);
      }
    }else{ // des is a node
      nodeSize2(ansIdx)+=nodeSize2(desIdx);
      if(node2cache(ansIdx-ntip)==-1){ //if the node is not in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        node2cache(ansIdx-ntip)=cacheIdx; // link the node to cache
        cacheFree(cacheIdx)=false; // mark cache as used
        nUse++;
      }
      int desCacheIdx=node2cache(desIdx-ntip);
      int ansCacheIdx=node2cache(ansIdx-ntip);
      for(int x=0;x<nbin;x++){ // update member of ans
        cacheMat(x,ansCacheIdx)|=cacheMat(x,desCacheIdx);
      }
      
      //desIdx search for similar nodes in edge1
      int nodeStop=startNode(desIdx-ntip,0); //the closest node in edge1
      int xorStop=ntip*2; //the closest distance. Initial value is the maximal distance
      for(int j=0;j<2;j++){
        int currentNode=startNode(desIdx-ntip,j); //search begin from two startNode
        int xorDis=0;
        int nLarger=0; //nLarger is the No. tips only in node of edge1
        while(currentNode>ntip&&nLarger<=nodeSize2(desIdx)*threshold){ //if currentNode does not reach root and currentNode is not much larger than desIdx
          IntegerVector tmpDiff=countDifferentBits(mat[nnode1-1-(currentNode-ntip)],
                                                   cacheMat.column(desCacheIdx));
          xorDis=tmpDiff(0);
          nLarger=tmpDiff(1);
          int tmpIdx=nnode1-1-(currentNode-ntip);
          if(outputDis1(tmpIdx)>=xorDis){ //update min distance of nodes in edge1 to nodes in edge2, >= means preferring nodes find later (larger nodes)
            outputDis1(tmpIdx)=xorDis;
            outputNode1(tmpIdx)=desIdx+1;
          }
          if(xorStop>=xorDis){ //update closest node in edge1 to desIdx, >= means preferring nodes find later (larger nodes)
            xorStop=xorDis;
            nodeStop=currentNode;
            outputDis2(desIdx-ntip)=xorDis;
            outputNode2(desIdx-ntip)=currentNode+1;
          }
          //move currentNode to upper and start new comparison
          currentNode=linkVector(currentNode);
        }
      }
      if(startNode(ansIdx-ntip,0)==-1){
        startNode(ansIdx-ntip,0)=nodeStop;
      }else{
        startNode(ansIdx-ntip,1)=nodeStop;
      }
      //search finished
      
      cacheFree(desCacheIdx)=true; //unlink des from cache
      nUse--;
      for(int x=0;x<nbin;x++){ // clean cache
        cacheMat(x,desCacheIdx)=0;
      }
    }
  }
  IntegerMatrix output1(nnode1-1,2);
  IntegerMatrix output2(nnode2,3);
  output1.column(0)=outputNode1;
  output1.column(1)=outputDis1;
  output2.column(0)=outputNode2;
  output2.column(1)=outputDis2; 
  for(int i=0;i<nnode2;i++){
    output2(i,2)=nodeSize2(i+ntip);
  }
  return Rcpp::List::create(
    Named("ref2qry")=output1,
    Named("qry2ref")=output2
  );
}

// [[Rcpp::export]]
std::vector<int> findMaxDes(IntegerMatrix edge,int root,IntegerVector nodeState){
  // for a given node (root), find largest descendants (nodeState>0)
  std::vector<int> output;
  std::vector<int> ances_stack;
  ances_stack.push_back(root);
  int i=0;
  //search for the root
  while(edge(i,0)!=root&&i<edge.nrow()){
    i++;
  }
  //search for descendants
  while(i<edge.nrow()&&ances_stack.size()>0){
    int ansIdx=edge(i,0);
    int desIdx=edge(i,1);
    if(ansIdx==ances_stack.back()){//if ansIdx match current ancestor
      if(nodeState(desIdx-1)>0){//if desIdx is one wanted node
        output.push_back(desIdx);
      }else{//if not, change the current ancestor
        ances_stack.push_back(desIdx);
      }
      i++;
    }else if(ansIdx<ances_stack.back()){// ansIdx returns to higher node
      ances_stack.pop_back();
    }else{// ansIdx is still under current ancestor
      i++;
    }
  }
  return output;
}

// [[Rcpp::export]]
IntegerMatrix qry2ref(IntegerVector nodes,IntegerVector xorDis,int ntip,int nnode){
  //transform data format from: each node in ref to:each node in qry
  //ntip: No. tips nnode No. nodes in qry
  
  IntegerVector outputDis(nnode);
  outputDis.fill(ntip*2);
  IntegerVector outputNode(nnode);
  outputNode.fill(-1);
  for(int i=0;i<xorDis.size();i++){
    int idx=nodes(i)-ntip-1;
    if(idx>=0){
      if(outputDis(idx)>xorDis(i)){
        outputDis(idx)=xorDis(i);
        outputNode(idx)=i+1+ntip;
      }
    }
  }
  IntegerMatrix output(nnode,2);
  output.column(0)=outputNode;
  output.column(1)=outputDis;
  return output;
}

// [[Rcpp::export]]
IntegerVector collapseNode(IntegerMatrix edge,IntegerVector nodeState,int ntip){
  IntegerMatrix output=clone(edge);
  IntegerVector linkVector=edge2vector(edge,ntip*2+1);
  for(int i=0;i<output.nrow();i++){
    int ansIdx=output(i,0)-1;
    int desIdx=output(i,1)-1;
    if(nodeState(desIdx)==0){
      output(i,0)=-1;
      output(i,1)=-1;
    }else{
      while(nodeState(ansIdx)==0){
        ansIdx=linkVector(ansIdx);
      }
      output(i,0)=ansIdx+1;
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> removeNestedNode(IntegerMatrix edge,IntegerVector nodes,int size){
  IntegerVector linkVector=edge2vector(edge,size);
  IntegerVector nodesVector(size);
  nodesVector.fill(0);
  for(int i=0;i<nodes.size();i++){
    int idx=nodes(i)-1;
    nodesVector(idx)=1;
  }
  std::vector<int> output;
  for(int i=0;i<nodes.size();i++){
    int idx=nodes(i)-1;
    int currentNode=linkVector(idx);
    while(currentNode!=-1){
      if(nodesVector(currentNode)==1){
        break;
      }
      currentNode=linkVector(currentNode);
    }
    if(currentNode==-1){
      output.push_back(idx+1);
    }
  }
  return output;
}

// [[Rcpp::export]]
std::vector<int> buildNode(IntegerMatrix edge,IntegerVector nodeSize,IntegerVector nodes,int ntip,float threshold=0.1){
  int maxIdx=nodeSize.size();
  IntegerVector linkVector=edge2vector(edge,maxIdx);
  IntegerVector inCount(maxIdx); //No. tips included in the group
  inCount.fill(0);
  IntegerVector passCount(maxIdx); //when search from nodes to root, No. times passed it.
  passCount.fill(0);
  int totalTips=0; //total tips from all nodes to join
  
  //prepare inCount and passCount
  for(int i=0;i<nodes.size();i++){
    int currentNode=nodes(i)-1;
    totalTips+=nodeSize(currentNode); //count total tips
    int inSize=nodeSize(currentNode);
    do{
      inCount(currentNode)+=inSize;
      passCount(currentNode)++;
      currentNode=linkVector(currentNode);
    } while (currentNode!=-1);
  }
  
  std::vector<int> nodesRun=Rcpp::as<std::vector<int>>(nodes); //candidate nodes to rebuild a clade 
  while(nodesRun.size()>1){
    std::vector<int> nodesMerge; //the ancestors input node attached to, these ancestors are direct decesdants of MRCA and can be joined together
    IntegerVector included(maxIdx); //whether the node is included in nodesMerge
    included.fill(0);
    IntegerVector nodesAns(nodesRun.size()); //for each input nodes,which node it is attached to
    nodesAns.fill(-1);
    for(unsigned int i=0;i<nodesRun.size();i++){
      int currentNode=nodesRun[i]-1;
      do{
        int tmp=linkVector(currentNode);
        if((unsigned int)(passCount(tmp))>=nodesRun.size()){ // tmp is the most recent common ancestors of of nodesRun (all nodes passed it)
          nodesAns(i)=currentNode; //ith input node is attached to the node found
          if(included(currentNode)==0){ //if the node is not included
            nodesMerge.push_back(currentNode);
            included(currentNode)=1;
          }
          break;
        }else{ //if not found, move to ancestor node
          currentNode=tmp;        
        }
      } while (currentNode!=-1);
    }
    //sort nodeMerge by the value nodeSize-inCount*2
    //when including node i, deltaDis=(nodeSize(i)-inCount(i))-inCount(i)
    //So sort nodes by deltaDis
    std::sort(
      nodesMerge.begin(),
      nodesMerge.end(),
      [&](int i, int j){
        return nodeSize(i)-inCount(i)*2<nodeSize(j)-inCount(j)*2;
      }
    );
    std::vector<int> output;
    int totalDelta=0; //delta distance when some nodes are included
    int subTotalTips=0; //total tips joined together
    int tmpAns=linkVector(nodesMerge[0]); //the MRCA of merged tips
    for(unsigned int i=0;i<nodesMerge.size();i++){
      int tmpNode=nodesMerge[i]; // the node index
      int tmpDelta=nodeSize(tmpNode)-inCount(tmpNode)*2; //when including the node, the delta of distance
      if(tmpDelta<=0){ //if delta<=0, including the node decrease the distance
        totalDelta+=tmpDelta;
        output.push_back(tmpNode+1);
        subTotalTips+=nodeSize(tmpNode);
      }else{ //if delta>0, including the node increase the distance
        break;
      }
    }
    //totalTips+totalDelta is the final distance
    if((totalTips+totalDelta)<=totalTips*threshold){ // find some nodes to rebuild a similar clade
      if(subTotalTips==nodeSize(tmpAns)){ //if all nodes joined rebuild their MRCA, return MRCA
        output.clear();
        output.push_back(tmpAns+1);
      }
      return(output);
    }else{
      int nextCommonAns=-1;
      for(unsigned int i=0;i<nodesMerge.size();i++){
        if(inCount(nodesMerge[i])>=totalTips*(1-threshold)){
          nextCommonAns=nodesMerge[i];
          break;
        }
      }
      if(nextCommonAns!=-1){ //prepare nodes for next round of search
        std::vector<int> newNodes;
        for(int i=0;i<nodesAns.size();i++){
          if(nodesAns[i]==nextCommonAns){
            newNodes.push_back(nodesRun[i]);
          }
        }
        nodesRun=newNodes;
      }else{ //cannot find a subclade for next round of search
        return(std::vector<int>());
      }
    }
  }
  // one node includes most tips from all nodes
  if(nodesRun.size()==1){
    int tmpNode=nodesRun[0]-1;
    if(totalTips*threshold>=totalTips-inCount(tmpNode)+nodeSize(tmpNode)-inCount(tmpNode)){
      return(nodesRun);
    }
  }
  return(std::vector<int>());
}

// [[Rcpp::export]]
Rcpp::List getDesIdx(IntegerMatrix edge,int ntip,int nnode){
  IntegerVector output(edge.nrow());
  IntegerVector nDes(nnode);
  nDes.fill(0);
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-ntip-1;
    nDes(ansIdx)++;
  }
  IntegerVector offset(nnode+1);
  offset.fill(0);
  for(int i=1;i<=nnode;++i){
    offset(i)=offset(i-1)+nDes(i-1);
  }
  nDes.fill(0);
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-ntip-1;
    int idx=nDes(ansIdx)+offset(ansIdx);
    output(idx)=edge(i,1);
    nDes(ansIdx)++;
  }
  return Rcpp::List::create(
    Named("idx")=output,
    Named("offset")=offset
  );
}

// [[Rcpp::export]]
IntegerVector getNodeDescendants(List desResult,int ntip,int node) {
  IntegerVector idx = desResult["idx"];
  IntegerVector offset = desResult["offset"];
  
  int ansIdx = node - ntip - 1;
  int start = offset[ansIdx];
  int numDescendants = offset[ansIdx+1]-offset[ansIdx];
  
  IntegerVector descendants(numDescendants);
  for(int i = 0; i < numDescendants; ++i) {
    descendants[i] = idx[start + i];
  }
  return descendants;
}

// [[Rcpp::export]]
Rcpp::List approximateCompare2(SEXP bitClade,int ntip,
                               IntegerMatrix edge1,
                               IntegerMatrix edge2,int nnode2,
                               int cacheSize=8192,float threshold=0.1){
  XPtr<BigMatrix> xpMat(bitClade);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int nnode1=xpMat->ncol()+1; //root node is not included in bitClade
  int nbin=xpMat->nrow();
  IntegerMatrix cacheMat(nbin,cacheSize);
  cacheMat.fill(0);
  LogicalVector cacheFree(cacheSize);
  cacheFree.fill(true);
  int nUse=0;
  IntegerVector node2cache(nnode2);  //ith element contains the cache idx of ith node 
  node2cache.fill(-1);
  int cacheIdx=0;
  
  IntegerVector outputDis1(nnode1-1);
  outputDis1.fill(ntip*2);
  IntegerVector outputNode1(nnode1-1);
  outputNode1.fill(0);
  IntegerVector outputDis2(nnode2);
  outputDis2.fill(ntip*2);
  IntegerVector outputNode2(nnode2);
  outputNode2.fill(0);
  
  //create search idx for edge1
  IntegerVector linkVector=edge2vector(edge1,ntip+nnode1);
  
  //create index to quick search of descendants in edge2  
  Rcpp::List desResult=getDesIdx(edge2,ntip,nnode2);
  
  //the node of edge1 where node of edge2 stop search
  IntegerVector nodeStop(ntip+nnode2);
  //ntip is the idx of root node in edge1
  nodeStop.fill(ntip);
  //the sizes of nodes in edge2
  IntegerVector nodeSize2(ntip+nnode2);
  nodeSize2.fill(0);
  
  for(int i=edge2.nrow()-1;i>=0;i--){
    int ansIdx=edge2(i,0)-1;
    int desIdx=edge2(i,1)-1;
    if(desIdx<ntip){ // des is a tip
      nodeSize2(desIdx)=1;
      nodeSize2(ansIdx)++;
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      int ansCacheIdx=node2cache(ansIdx-ntip); // find ans in cache
      if(ansCacheIdx==-1){ //if ans is not in cache,store it in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        ansCacheIdx=cacheIdx;
        node2cache(ansIdx-ntip)=cacheIdx; // link node to cache
        cacheFree(cacheIdx)=false; //mark cache as used
        nUse++;
      }
      cacheMat(x,ansCacheIdx)|=(1<<y); // update the member of ans
      // when ansIdx start searching, it can begin from the ancestor of desIdx in edge1
      nodeStop(desIdx)=linkVector(desIdx);
    }else{ // des is a node
      nodeSize2(ansIdx)+=nodeSize2(desIdx);
      if(node2cache(ansIdx-ntip)==-1){ //if the node is not in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        node2cache(ansIdx-ntip)=cacheIdx; // link the node to cache
        cacheFree(cacheIdx)=false; // mark cache as used
        nUse++;
      }
      int desCacheIdx=node2cache(desIdx-ntip);
      int ansCacheIdx=node2cache(ansIdx-ntip);
      for(int x=0;x<nbin;x++){ // update member of ans
        cacheMat(x,ansCacheIdx)|=cacheMat(x,desCacheIdx);
      }
      
      //desIdx search for similar nodes in edge1
      IntegerVector startNode=getNodeDescendants(desResult,ntip,desIdx+1);
      for(int j=0;j<startNode.size();++j){
        startNode(j)=nodeStop(startNode(j)-1);
      }
      nodeStop(desIdx)=startNode(0);
      int xorStop=ntip*2; //the closest distance. Initial value is the maximal distance
      for(int j=0;j<startNode.size();++j){
        int currentNode=startNode(j); //search begin from startNode
        int xorDis=0;
        int nLarger=0; //nLarger is the No. tips only in node of edge1
        while(currentNode>ntip&&nLarger<=nodeSize2(desIdx)*threshold){ //if currentNode does not reach root and currentNode is not much larger than desIdx
          IntegerVector tmpDiff=countDifferentBits(mat[nnode1-1-(currentNode-ntip)],
                                                   cacheMat.column(desCacheIdx));
          xorDis=tmpDiff(0);
          nLarger=tmpDiff(1);
          int tmpIdx=nnode1-1-(currentNode-ntip);
          if(outputDis1(tmpIdx)>=xorDis){ //update min distance of nodes in edge1 to nodes in edge2, >= means preferring nodes find later (larger nodes)
            outputDis1(tmpIdx)=xorDis;
            outputNode1(tmpIdx)=desIdx+1;
          }
          if(xorStop>=xorDis){ //update closest node in edge1 to desIdx, >= means preferring nodes find later (larger nodes)
            xorStop=xorDis;
            nodeStop(desIdx)=currentNode;
            outputDis2(desIdx-ntip)=xorDis;
            outputNode2(desIdx-ntip)=currentNode+1;
          }
          //move currentNode to upper and start new comparison
          currentNode=linkVector(currentNode);
        }
      }
      //search finished
      
      cacheFree(desCacheIdx)=true; //unlink des from cache
      nUse--;
      for(int x=0;x<nbin;x++){ // clean cache
        cacheMat(x,desCacheIdx)=0;
      }
    }
  }
  IntegerMatrix output1(nnode1,2);
  output1(0,0)=ntip+1;
  output1(0,1)=0;
  int k=nnode1-2;
  for(int i=1;i<nnode1;i++){
    output1(i,0)=outputNode1(k);
    output1(i,1)=outputDis1(k);
    k--;
  }  
  IntegerMatrix output2(nnode2,3);
  output2.column(0)=outputNode2;
  output2.column(1)=outputDis2; 
  for(int i=0;i<nnode2;i++){
    output2(i,2)=nodeSize2(i+ntip);
  }
  output2(0,0)=ntip+1;
  output2(0,1)=0;
  return Rcpp::List::create(
    Named("ref2qry")=output1,
    Named("qry2ref")=output2
  );
}

IntegerVector countDifferentBits_AVX2(int* xx, IntegerVector yy) {
  // count different bits (count(0)) and unique bits in xx (count(1))
  IntegerVector count(2);
  count.fill(0);
  
  size_t n = yy.size();
  size_t i = 0;
  
  __m256i sum_vec0 = _mm256_setzero_si256();
  __m256i sum_vec1 = _mm256_setzero_si256();
  
  for (; i + 8 <= n; i += 8) {
    __m256i va = _mm256_loadu_si256((const __m256i*)&xx[i]);
    __m256i vb = _mm256_loadu_si256((const __m256i*)&yy[i]); 
    
    __m256i vxor = _mm256_xor_si256(va, vb);
    __m256i pop0 = popcount_swar_avx2(vxor);
    sum_vec0 = _mm256_add_epi32(sum_vec0, pop0);
    
    __m256i vand = _mm256_andnot_si256(vb, va);  // va & (~vb)
    __m256i pop1 = popcount_swar_avx2(vand);
    sum_vec1 = _mm256_add_epi32(sum_vec1, pop1);
  }
  
  int sum_array0[8], sum_array1[8];
  _mm256_storeu_si256((__m256i*)sum_array0, sum_vec0);
  _mm256_storeu_si256((__m256i*)sum_array1, sum_vec1);
  
  int total_count0 = 0, total_count1 = 0;
  for (int j = 0; j < 8; j++) {
    total_count0 += sum_array0[j];
    total_count1 += sum_array1[j];
  }
  
  for (; i < n; ++i) {
    total_count0 += popcount_scalar(xx[i] ^ yy[i]);
    total_count1 += popcount_scalar(xx[i] & (~yy[i]));
  }
  
  count(0) = total_count0;
  count(1) = total_count1;
  
  return count;
}

// [[Rcpp::export]]
Rcpp::List approximateCompare2_AVX2(SEXP bitClade,int ntip,
                                    IntegerMatrix edge1,
                                    IntegerMatrix edge2,int nnode2,
                                    int cacheSize=8192,float threshold=0.1){
  XPtr<BigMatrix> xpMat(bitClade);
  MatrixAccessor<int> mat=MatrixAccessor<int>(*xpMat);
  
  int nnode1=xpMat->ncol()+1; //root node is not included in bitClade
  int nbin=xpMat->nrow();
  IntegerMatrix cacheMat(nbin,cacheSize);
  cacheMat.fill(0);
  LogicalVector cacheFree(cacheSize);
  cacheFree.fill(true);
  int nUse=0;
  IntegerVector node2cache(nnode2);  //ith element contains the cache idx of ith node 
  node2cache.fill(-1);
  int cacheIdx=0;
  
  IntegerVector outputDis1(nnode1-1);
  outputDis1.fill(ntip*2);
  IntegerVector outputNode1(nnode1-1);
  outputNode1.fill(0);
  IntegerVector outputDis2(nnode2);
  outputDis2.fill(ntip*2);
  IntegerVector outputNode2(nnode2);
  outputNode2.fill(0);
  
  //create search idx for edge1
  IntegerVector linkVector=edge2vector(edge1,ntip+nnode1);
  
  //create index to quick search of descendants in edge2  
  Rcpp::List desResult=getDesIdx(edge2,ntip,nnode2);
  
  //the node of edge1 where node of edge2 stop search
  IntegerVector nodeStop(ntip+nnode2);
  //ntip is the idx of root node in edge1
  nodeStop.fill(ntip);
  //the sizes of nodes in edge2
  IntegerVector nodeSize2(ntip+nnode2);
  nodeSize2.fill(0);
  
  for(int i=edge2.nrow()-1;i>=0;i--){
    int ansIdx=edge2(i,0)-1;
    int desIdx=edge2(i,1)-1;
    if(desIdx<ntip){ // des is a tip
      nodeSize2(desIdx)=1;
      nodeSize2(ansIdx)++;
      int x=desIdx>>5; // divided by 32
      int y=desIdx&31; // mod by 32
      int ansCacheIdx=node2cache(ansIdx-ntip); // find ans in cache
      if(ansCacheIdx==-1){ //if ans is not in cache,store it in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        ansCacheIdx=cacheIdx;
        node2cache(ansIdx-ntip)=cacheIdx; // link node to cache
        cacheFree(cacheIdx)=false; //mark cache as used
        nUse++;
      }
      cacheMat(x,ansCacheIdx)|=(1<<y); // update the member of ans
      // when ansIdx start searching, it can begin from the ancestor of desIdx in edge1
      nodeStop(desIdx)=linkVector(desIdx);
    }else{ // des is a node
      nodeSize2(ansIdx)+=nodeSize2(desIdx);
      if(node2cache(ansIdx-ntip)==-1){ //if the node is not in cache
        if(nUse==cacheSize){
          Rprintf("Run out of cache!\n");
          return 0;
        }
        while(!cacheFree(cacheIdx)){ //find a free cache
          cacheIdx++;
          if(cacheIdx>=cacheSize){
            cacheIdx=0;
          }
        }
        node2cache(ansIdx-ntip)=cacheIdx; // link the node to cache
        cacheFree(cacheIdx)=false; // mark cache as used
        nUse++;
      }
      int desCacheIdx=node2cache(desIdx-ntip);
      int ansCacheIdx=node2cache(ansIdx-ntip);
      for(int x=0;x<nbin;x++){ // update member of ans
        cacheMat(x,ansCacheIdx)|=cacheMat(x,desCacheIdx);
      }
      
      //desIdx search for similar nodes in edge1
      IntegerVector startNode=getNodeDescendants(desResult,ntip,desIdx+1);
      for(int j=0;j<startNode.size();++j){
        startNode(j)=nodeStop(startNode(j)-1);
      }
      nodeStop(desIdx)=startNode(0);
      int xorStop=ntip*2; //the closest distance. Initial value is the maximal distance
      for(int j=0;j<startNode.size();++j){
        int currentNode=startNode(j); //search begin from startNode
        int xorDis=0;
        int nLarger=0; //nLarger is the No. tips only in node of edge1
        while(currentNode>ntip&&nLarger<=nodeSize2(desIdx)*threshold){ //if currentNode does not reach root and currentNode is not much larger than desIdx
          IntegerVector tmpDiff=countDifferentBits_AVX2(mat[nnode1-1-(currentNode-ntip)],
                                                        cacheMat.column(desCacheIdx));
          xorDis=tmpDiff(0);
          nLarger=tmpDiff(1);
          int tmpIdx=nnode1-1-(currentNode-ntip);
          if(outputDis1(tmpIdx)>=xorDis){ //update min distance of nodes in edge1 to nodes in edge2, >= means preferring nodes find later (larger nodes)
            outputDis1(tmpIdx)=xorDis;
            outputNode1(tmpIdx)=desIdx+1;
          }
          if(xorStop>=xorDis){ //update closest node in edge1 to desIdx, >= means preferring nodes find later (larger nodes)
            xorStop=xorDis;
            nodeStop(desIdx)=currentNode;
            outputDis2(desIdx-ntip)=xorDis;
            outputNode2(desIdx-ntip)=currentNode+1;
          }
          //move currentNode to upper and start new comparison
          currentNode=linkVector(currentNode);
        }
      }
      //search finished
      
      cacheFree(desCacheIdx)=true; //unlink des from cache
      nUse--;
      for(int x=0;x<nbin;x++){ // clean cache
        cacheMat(x,desCacheIdx)=0;
      }
    }
  }
  IntegerMatrix output1(nnode1,2);
  output1(0,0)=ntip+1;
  output1(0,1)=0;
  int k=nnode1-2;
  for(int i=1;i<nnode1;i++){
    output1(i,0)=outputNode1(k);
    output1(i,1)=outputDis1(k);
    k--;
  }  
  IntegerMatrix output2(nnode2,3);
  output2.column(0)=outputNode2;
  output2.column(1)=outputDis2; 
  for(int i=0;i<nnode2;i++){
    output2(i,2)=nodeSize2(i+ntip);
  }
  output2(0,0)=ntip+1;
  output2(0,1)=0;
  return Rcpp::List::create(
    Named("ref2qry")=output1,
    Named("qry2ref")=output2
  );
}

// [[Rcpp::export]]
std::vector<int> findClone(IntegerMatrix& edge,int ntip,int nnode,NumericVector& bs,double threshold){
  IntegerVector nodeState(nnode);
  nodeState.fill(1);
  nodeState[0]=0;
  std::vector<int> output;
  for(int i=0;i<edge.nrow();++i){
    int ansIdx=edge(i,0)-ntip-1;
    int desIdx=edge(i,1)-ntip-1;
    if(nodeState[ansIdx]==0&&desIdx>=0){
      if(bs[desIdx]>=threshold){
        output.push_back(desIdx+1+ntip);
      }else{
        nodeState[desIdx]=0;
      }
    }
  }
  return output;
}
