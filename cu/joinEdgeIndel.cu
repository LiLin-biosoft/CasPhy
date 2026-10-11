#include "cuda_runtime.h"
#include <cuda_fp16.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define BLOCKSIZE 128
#define MAXCOL 1000

#define CUDACALL(call) {                                    \
	int nTry=0;	\
	while(nTry<10){ \
		nTry++;                                                    \
		cudaError_t err = call;                                         \
		if(err==cudaSuccess){ \
			break; \
		}else{ \
			printf("CUDA error at %s:%d code=%d(%s) \"%s.\" Try %d times.\n",          \
				__FILE__, __LINE__, err, cudaGetErrorString(err), #call,nTry); \
			sleep(60); \
			if(nTry>=10){ \
				exit(1); \
			} \
		}                                                               \
	} \
}

//joinEdgeIndelV2

typedef unsigned int u32;

__constant__ u32 reference[MAXCOL];

__global__ void fillFloat(half* const data,u32 const size,half const value){
	u32 tid=blockIdx.x*blockDim.x+threadIdx.x;
	u32 stride=blockDim.x*gridDim.x;
    for(u32 i=tid;i<size;i+=stride){
		data[i]=value;
    }
}

__global__ void fillFloatCol(half* const data,u32 const size,half const value){
	u32 tid=blockIdx.x*blockDim.x+threadIdx.x;
	u32 stride=blockDim.x*gridDim.x;
    for(u32 i=tid;i<size;i+=stride){
		data[i*size]=value;
    }
}

__global__ void FP16toFP32Kernel(float* const output,half* input,u32 const size){
	u32 tid=blockIdx.x*blockDim.x+threadIdx.x;
	u32 stride=blockDim.x*gridDim.x;
    for(u32 i=tid;i<size;i+=stride){
		output[i]=__half2float(input[i]);
    }
}

__global__ void FP32toFP16Kernel(half* const output,float* input,u32 const size){
	u32 tid=blockIdx.x*blockDim.x+threadIdx.x;
	u32 stride=blockDim.x*gridDim.x;
    for(u32 i=tid;i<size;i+=stride){
		output[i]=__float2half(input[i]);
    }
}

__global__ void calculateSharedKernel(u32* const cacheMat,u32 const ncol, 
							u32* const rowIndices,u32 nIndices,
							half* const cacheShared,half* const weight, u32 const selfIdx){
	const u32 tid=threadIdx.x+blockIdx.x*blockDim.x;
	if(tid<nIndices){
		const u32 rowIdx=rowIndices[tid];
		half commonCount=0;
		for(u32 i=0;i<ncol;i++){
			u32 tmp=reference[i];
			if(tmp==cacheMat[rowIdx*ncol+i]&&tmp>1){
				commonCount+=weight[(tmp-2)*ncol+i]; //row: state,col: site
			}
		}
		if(rowIdx==selfIdx){
			commonCount=-1;
		}
		cacheShared[rowIdx]=commonCount;
	}
}

__global__ void updateAncestralKernel(u32* const ve1,u32* const ve2,const u32 size){
	u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
	u32 stride = blockDim.x * gridDim.x;
	for (u32 i = tid; i < size; i += stride) {
		u32 tmp1=ve1[i];
		u32 tmp2=ve2[i];
		if(tmp1==0){
			ve1[i]=tmp2;
		}else if(tmp2>=1&&tmp1!=tmp2){
			ve1[i]=1;
		}
	}
}

__global__ void findMaxFloatKernel(const half* array,u32* const arrayIndices,u32 size,
							half* max_values,u32* max_indices){
    __shared__ half threadMax[BLOCKSIZE];
    __shared__ u32 threadMaxIdx[BLOCKSIZE];
    
    u32 tid=threadIdx.x;
    u32 global_id=blockIdx.x*blockDim.x+threadIdx.x;
    u32 stride=blockDim.x*gridDim.x;
    
    half local_max=-1;
    u32 local_index=0;
    
    // each thread process several elements
    for(u32 i=global_id;i<size;i+=stride){
    	half value=array[i];
    	u32 idx=arrayIndices[i];
       if(value>local_max||value==local_max&&idx<local_index){
            local_max=value;
            local_index=idx;
        }
    }
    
    // each thread save local max to shared memory
    threadMax[tid]=local_max;
    threadMaxIdx[tid]=local_index;
    __syncthreads();
    
    // find local max in each block
    for(u32 s=blockDim.x/2;s>0;s>>=1){
        if(tid<s){
            if(threadMax[tid+s]>threadMax[tid]||(threadMax[tid+s]==threadMax[tid]&&threadMaxIdx[tid+s]<threadMaxIdx[tid])) {
				threadMax[tid]=threadMax[tid+s];
                threadMaxIdx[tid]=threadMaxIdx[tid+s];
            }
        }
        __syncthreads();
    }
    
    // save block local max to output
    if(tid==0){
        max_values[blockIdx.x]=threadMax[0];
        max_indices[blockIdx.x]=threadMaxIdx[0];
    }
}

__global__ void findValueKernel(u32* const edge,u32 nedge,
							u32 value, u32* idx){
    __shared__ u32 threadFound[BLOCKSIZE];
    
    u32 tid=threadIdx.x;
    u32 stride=blockDim.x;
    
    u32 lineIdx=nedge+1;
    
    // each thread process several elements
    for(u32 i=tid;i<nedge;i+=stride){
		if(edge[i*2]==value&&i<=lineIdx){
			lineIdx=i;
		}
    }
    
    // each thread save local minIdx to shared memory
    threadFound[tid]=lineIdx;
    __syncthreads();
    
    // find local minIdx in each block
    for(u32 s=blockDim.x/2;s>0;s>>=1){
        if(tid<s){
            if(threadFound[tid]>threadFound[tid+s]) {
				threadFound[tid]=threadFound[tid+s];
            }
        }
        __syncthreads();
    }
    u32 globalFound=threadFound[0];
    
    // save global minIdx to output
    if(tid==0){
        idx[0]=globalFound;
	if(globalFound<nedge){
        	edge[globalFound*2]=(nedge+1);
	}
    }
}

__global__ void findValueKernel2(u32* const edge,u32 nedge,
							u32 value, u32* idx){
    u32 tid=threadIdx.x;
    u32 stride=blockDim.x;
    for(u32 i=tid;i<nedge;i+=stride){
		if(edge[i*2+1]==value){
			idx[0]=i;
		}
    }
}

__host__ void calculateSharedPacked(
				u32* rowIndices,u32 nIndices,
				u32* cacheMat,
				half* cacheShared_gpu,
				u32 ntip,u32 ncol,half* weight_gpu,
				u32 lineIdx
){
	u32* rowIndices_gpu;
	cudaMalloc((void**)&rowIndices_gpu,sizeof(u32)*ntip);

	// copy the indices of rows for calculation to GPU
	u32 BlockN=(nIndices-1)/BLOCKSIZE+1;
	cudaMemcpy(rowIndices_gpu,rowIndices,sizeof(u32)*nIndices,cudaMemcpyHostToDevice);
	// copy reference row to constant memory
	cudaMemcpyToSymbol(reference,cacheMat+lineIdx*ncol,sizeof(u32)*ncol);
	// calcuate shared
	fillFloat<<<1,BLOCKSIZE>>>(cacheShared_gpu,ntip,-1);
	calculateSharedKernel<<<BlockN,BLOCKSIZE>>>(
		cacheMat,ncol,
		rowIndices_gpu,nIndices,
		cacheShared_gpu,weight_gpu,lineIdx
	);

	cudaFree(rowIndices_gpu);
}

__host__ void findHalfMax(half* array,u32* idx,u32 nValue,half* maxValue,u32* maxIdx){
	u32 BlockN=(nValue-1)/BLOCKSIZE+1;
	half* maxValues_gpu;
	cudaMalloc((void**)&maxValues_gpu,sizeof(half)*BlockN);
	u32* maxIndices_gpu;
	cudaMalloc((void**)&maxIndices_gpu,sizeof(u32)*BlockN);
	findMaxFloatKernel<<<BlockN,BLOCKSIZE>>>(
		array,idx,nValue,
		maxValues_gpu,maxIndices_gpu
	);
	findMaxFloatKernel<<<1,BLOCKSIZE>>>(
		maxValues_gpu,maxIndices_gpu,BlockN,
		maxValue,maxIdx
	);
	cudaFree(maxValues_gpu);
	cudaFree(maxIndices_gpu);
}

//each node calculate shared mutations with previous nodes
//calculate cacheShared when using

int main(int argc, char *argv[]){
	setbuf(stdout, NULL); 
	clock_t start = clock();

	int gpuCount;
	CUDACALL(cudaGetDeviceCount(&gpuCount));
	srand(time(nullptr));
	int runDevice=rand()%gpuCount;
	CUDACALL(cudaSetDevice(runDevice));

	u32 minSize=atoi(argv[4]);
	float minRatio=atof(argv[5]);

	u32 nState=atoi(argv[3]);	
	FILE *file;
	char input_file[256];
	strcpy(input_file,argv[1]);
	strcat(input_file,".tbmat"); //read mat,each row is one tip
    u32 dim[2];
    file = fopen(input_file, "rb");
    if (file == NULL) {
    	printf("No tbmat input.\n");
    	return 0;
  	} 
    fread(dim, sizeof(u32), 2, file);
    const u32 ntip=dim[0],ncol=dim[1];
    const u32 blockN_ntip=(ntip-1)/BLOCKSIZE+1;
    const u32 nnode=ntip*2-1;
    const u32 nedge=ntip*2-2;
    const half negOne=-1;

	u32* cacheMatHost;
	cudaHostAlloc((void**)&cacheMatHost,sizeof(u32)*ncol*ntip,cudaHostAllocDefault);
	u32* cacheMat;
	cudaMalloc((void**)&cacheMat,sizeof(u32)*ntip*ncol);

	fread(cacheMatHost, sizeof(u32), ncol*ntip, file);
	cudaMemcpy(cacheMat,cacheMatHost, sizeof(u32)*ncol*ntip, cudaMemcpyHostToDevice);
	cudaFreeHost(cacheMatHost);
	fclose(file);
	
    file=fopen(argv[2], "rb"); //read weight, double vector
    if (file == NULL) {
    	printf("No weight input.\n");
    	return 0;
  	} 
  	
	// transform weight from double to half
  	double* weightTmp=(double*)(malloc(sizeof(double)*ncol*nState));
	fread(weightTmp, sizeof(double),ncol*nState, file);   
	fclose(file);	
	half* weight=(half*)(malloc(sizeof(half)*ncol*nState));
	for(u32 i=0;i<ncol*nState;i++){
		weight[i]=(half)(weightTmp[i]);
	}
	free(weightTmp);

	half* weight_gpu;
	cudaMalloc((void**)&weight_gpu,sizeof(half)*ncol*nState);
	cudaMemcpy(weight_gpu,weight,sizeof(half)*ncol*nState,cudaMemcpyHostToDevice);
	free(weight); 

	float* cacheShared=(float*)(malloc(sizeof(float)*ntip));
	half* cacheShared_gpuFP16;
	cudaMalloc((void**)&cacheShared_gpuFP16,sizeof(half)*ntip);
	float* cacheShared_gpuFP32;
	cudaMalloc((void**)&cacheShared_gpuFP32,sizeof(float)*ntip);
	half* localMaxFP16;
	cudaMalloc((void**)&localMaxFP16,sizeof(half)*ntip);
	float* localMaxFP32;
	cudaMalloc((void**)&localMaxFP32,sizeof(float)*ntip);
	half* localMaxHostFP16;
	cudaHostAlloc((void**)&localMaxHostFP16,sizeof(half)*ntip,cudaHostAllocDefault);
	float* localMaxHostFP32;
	cudaHostAlloc((void**)&localMaxHostFP32,sizeof(float)*ntip,cudaHostAllocDefault);
	u32* localMaxIdx=(u32*)(malloc(sizeof(u32)*ntip));

	u32* nodeState;
	cudaHostAlloc((void**)&nodeState,sizeof(u32)*ncol*nnode,cudaHostAllocDefault);
	u32* node2cache=(u32*)(malloc(sizeof(u32)*nnode));
	u32* cache2node;
	cudaMalloc((void**)&cache2node,sizeof(u32)*ntip);
	u32* validCacheIndices=(u32*)(malloc(sizeof(u32)*ntip));
	u32 nValidCache=ntip;
	// initialize links from node to cache
	node2cache[0]=0;
	validCacheIndices[0]=0;
	for(u32 i=1;i<ntip;i++){
		node2cache[i]=i;
		node2cache[i-1+ntip]=ntip; //=ntip means invalid
		validCacheIndices[i]=i;
	}
	// initialize links from cache to node
	cudaMemcpy(cache2node,node2cache,sizeof(u32)*ntip,cudaMemcpyHostToDevice);
	
	half* maxValuesGlobal_gpu;
	cudaMalloc((void**)&maxValuesGlobal_gpu,sizeof(half));
	u32* maxIndicesGlobal_gpu;
	cudaMalloc((void**)&maxIndicesGlobal_gpu,sizeof(u32));
	u32* nodeSize=(u32*)(malloc(sizeof(u32)*ntip*2));
	strcpy(input_file,argv[1]);
	strcat(input_file,".tipSize"); //read mat,each row is one tip
    file = fopen(input_file, "rb");
    if (file == NULL) {
    	printf("No tipSize input.\n");
    	return 0;
  	} 
	fread(dim, sizeof(u32), 2, file);
	fread(nodeSize,sizeof(u32),ntip,file);
	fclose(file);

	for(u32 i=0;i<ntip-1;i++){
		u32 nIndices=ntip-i-1;
		calculateSharedPacked(
				node2cache+i+1,nIndices,
				cacheMat,
				cacheShared_gpuFP16,
				ntip,ncol,weight_gpu,
				i
		);
		findHalfMax(cacheShared_gpuFP16+i+1,cache2node+i+1,nIndices,maxValuesGlobal_gpu,maxIndicesGlobal_gpu);
		cudaMemcpy(localMaxIdx+i,maxIndicesGlobal_gpu,sizeof(u32),cudaMemcpyDeviceToHost);
		cudaMemcpy(localMaxFP16+i,maxValuesGlobal_gpu,sizeof(half),cudaMemcpyDeviceToDevice);
	}
	// write the first row local max to GPU (invalid)
	cudaMemcpy(localMaxFP16+ntip-1,&negOne,sizeof(half),cudaMemcpyHostToDevice);

	// initialize edge and node index
	u32* edge=(u32*)(malloc(sizeof(u32)*nedge*2));
	u32 nodeCount=ntip;
	// start join
	
	while(nodeCount<nnode){
		u32 joinIdx1=0,joinIdx2=0;
		u32 cacheIdx1=0,cacheIdx2=0;
		do{  //find global maximum
			findHalfMax(localMaxFP16,cache2node,ntip,maxValuesGlobal_gpu,maxIndicesGlobal_gpu);
			cudaMemcpy(&joinIdx1,maxIndicesGlobal_gpu,sizeof(u32),cudaMemcpyDeviceToHost);
			
			cacheIdx1=node2cache[joinIdx1];
			joinIdx2=localMaxIdx[cacheIdx1];
			cacheIdx2=node2cache[joinIdx2];
		
			if(cacheIdx2>=ntip){ // if joinIdx2 is not in cache (has been used), update joinIdx1's local maximum
				// global max found but joinIdx2 is invalid
				calculateSharedPacked(
					validCacheIndices,nValidCache,
					cacheMat,
					cacheShared_gpuFP16,
					ntip,ncol,weight_gpu,
					cacheIdx1
				);
				findHalfMax(cacheShared_gpuFP16,cache2node,ntip,maxValuesGlobal_gpu,maxIndicesGlobal_gpu);
				cudaMemcpy(localMaxIdx+cacheIdx1,maxIndicesGlobal_gpu,sizeof(u32),cudaMemcpyDeviceToHost);
				cudaMemcpy(localMaxFP16+cacheIdx1,maxValuesGlobal_gpu,sizeof(half),cudaMemcpyDeviceToDevice);
			}
		}while(cacheIdx2>=ntip); // Loop until joinIdx2 is not used
		
		// write state of join node to edge,then update ancestral
		cacheIdx2=node2cache[joinIdx2];
		cudaMemcpy(nodeState+(nodeCount-ntip)*2*ncol,cacheMat+ncol*cacheIdx1,ncol*sizeof(u32),cudaMemcpyDeviceToHost);
		cudaMemcpy(nodeState+((nodeCount-ntip)*2+1)*ncol,cacheMat+ncol*cacheIdx2,ncol*sizeof(u32),cudaMemcpyDeviceToHost);
		
		if(nodeSize[joinIdx1]<minSize&&nodeSize[joinIdx1]*minRatio<nodeSize[joinIdx2]){
			cudaMemcpy(cacheMat+ncol*cacheIdx1,cacheMat+ncol*cacheIdx2,ncol*sizeof(u32),cudaMemcpyDeviceToDevice);
		}else if(nodeSize[joinIdx2]<minSize&&nodeSize[joinIdx2]*minRatio<nodeSize[joinIdx1]){
			;
		}else{
			updateAncestralKernel<<<2,BLOCKSIZE>>>(cacheMat+ncol*cacheIdx1,cacheMat+ncol*cacheIdx2,ncol);
		}
		nodeSize[nodeCount]=nodeSize[joinIdx1]+nodeSize[joinIdx2];
		
		// write new branches to edge
		edge[(nodeCount-ntip)*4]=nodeCount;
		edge[(nodeCount-ntip)*4+1]=joinIdx1;
		edge[(nodeCount-ntip)*4+2]=nodeCount;
		edge[(nodeCount-ntip)*4+3]=joinIdx2;	
	
		// update index vector
		node2cache[joinIdx1]=ntip; //mark joinIdx1 as used
		node2cache[joinIdx2]=ntip; //mark joinIdx2 as used
		node2cache[nodeCount]=cacheIdx1; // link new node to cacheIdx1
		cudaMemcpy(cache2node+cacheIdx1,&nodeCount,sizeof(u32),cudaMemcpyHostToDevice);//link cacheIdx1 to new node
		cudaMemcpy(cache2node+cacheIdx2,&nnode,sizeof(u32),cudaMemcpyHostToDevice);//mark cacheIdx2 as invalid
		
		// update validCacheIndices
		nValidCache--;
		u32 ptr;
		for(ptr=0;validCacheIndices[ptr]<cacheIdx2;ptr++){
			;
		}
		for(;ptr<nValidCache;ptr++){
			validCacheIndices[ptr]=validCacheIndices[ptr+1];
		}
		
		calculateSharedPacked(
				validCacheIndices,nValidCache,
				cacheMat,
				cacheShared_gpuFP16,
				ntip,ncol,weight_gpu,
				cacheIdx1
		);

		FP16toFP32Kernel<<<blockN_ntip,BLOCKSIZE>>>(localMaxFP32,localMaxFP16,ntip);
		FP16toFP32Kernel<<<blockN_ntip,BLOCKSIZE>>>(cacheShared_gpuFP32,cacheShared_gpuFP16,ntip);
		cudaMemcpy(cacheShared,cacheShared_gpuFP32,sizeof(float)*ntip,cudaMemcpyDeviceToHost);
		cudaMemcpy(localMaxHostFP16,localMaxFP16,sizeof(half)*ntip,cudaMemcpyDeviceToHost);
		cudaMemcpy(localMaxHostFP32,localMaxFP32,sizeof(float)*ntip,cudaMemcpyDeviceToHost);

		for(u32 i=0;i<nValidCache;i++){
			u32 iTmp=validCacheIndices[i];
			float newShared=cacheShared[iTmp];
			if(cacheIdx1!=iTmp&&newShared>localMaxHostFP32[iTmp]){
				localMaxHostFP16[iTmp]=__float2half(newShared);
				localMaxIdx[iTmp]=nodeCount;
			}
		}

		localMaxHostFP16[cacheIdx1]=negOne;
		localMaxHostFP16[cacheIdx2]=negOne;
		cudaMemcpy(localMaxFP16,localMaxHostFP16,sizeof(half)*ntip,cudaMemcpyHostToDevice);


		nodeCount++;
	}
	cudaMemcpy(nodeState+(nnode-1)*ncol,cacheMat+ncol*node2cache[nodeCount-1],sizeof(u32)*ncol,cudaMemcpyDeviceToHost);
	free(nodeSize);
	cudaFree(cacheShared_gpuFP16);
	cudaFree(cacheShared_gpuFP32);
	cudaFree(cacheMat);
	cudaFree(weight_gpu);
	cudaFree(localMaxFP16);
	cudaFreeHost(localMaxHostFP16);
	cudaFree(localMaxFP32);
	cudaFreeHost(localMaxHostFP32);
	free(localMaxIdx);
	free(node2cache);
	cudaFree(cache2node);
	free(validCacheIndices);
	cudaFree(maxValuesGlobal_gpu);
	cudaFree(maxIndicesGlobal_gpu);

//////////////////////////////////////////

	u32* edge_gpu;
	cudaMalloc((void**)&edge_gpu,sizeof(u32)*nedge*2);
	cudaMemcpy(edge_gpu,edge, sizeof(u32)*nedge*2, cudaMemcpyHostToDevice);
	
	const u32 root=nedge; //ntip=(nedge+2)/2
	u32* nodeRelabel=(u32*)(malloc(sizeof(u32)*nnode));
	for(u32 i=0;i<ntip;i++){
		nodeRelabel[i]=i;
	}
	for(u32 i=ntip;i<nnode;i++){
		nodeRelabel[i]=nnode;
	}
	nodeCount=ntip;
	u32 idx=0;
	u32* idx_gpu;
	cudaMalloc((void**)&idx_gpu,sizeof(u32));	
	u32* ances_stack=(u32*)(malloc(sizeof(u32)*(nedge+1)));
	ances_stack[0]=root;
	int stack_pointer=0;

	FILE *outputEdge;
	char output_file[256];
	strcpy(output_file,argv[1]);
	strcat(output_file,".edge2");
	outputEdge = fopen(output_file, "wb");
	dim[0]=nedge;
	dim[1]=2;
	fwrite(dim, sizeof(u32), 2, outputEdge);

	FILE *outputNode;
	strcpy(output_file,argv[1]);
	strcat(output_file,".nodeState2");
	outputNode = fopen(output_file, "wb");
	dim[0]=ntip-1;
	dim[1]=ncol;
	fwrite(dim, sizeof(u32), 2, outputNode);
	fwrite(nodeState+(nedge)*ncol,sizeof(u32),ncol,outputNode);
	
	while(stack_pointer>(-1)){
		findValueKernel<<<1,BLOCKSIZE>>>(edge_gpu,nedge,ances_stack[stack_pointer],idx_gpu); 
		cudaMemcpy(&idx,idx_gpu,sizeof(u32),cudaMemcpyDeviceToHost);
		if(idx>=nedge){
			stack_pointer--;
		}else{
			stack_pointer++;
			ances_stack[stack_pointer]=edge[idx*2+1];
			u32 tmp[2];
			tmp[0]=edge[idx*2];tmp[1]=edge[idx*2+1];
			
			if(nodeRelabel[tmp[0]]==nnode){
				nodeRelabel[tmp[0]]=nodeCount;
				tmp[0]=nodeCount;
				nodeCount++;
			}else{
				tmp[0]=nodeRelabel[tmp[0]];
			}
			if(nodeRelabel[tmp[1]]==nnode){
				findValueKernel2<<<1,BLOCKSIZE>>>(edge_gpu,nedge,tmp[1],idx_gpu); 
				cudaMemcpy(&idx,idx_gpu,sizeof(u32),cudaMemcpyDeviceToHost);
//				printf("%u,%u\n",tmp[1],idx+1);
				fwrite(nodeState+idx*ncol,sizeof(u32),ncol,outputNode);
				
				nodeRelabel[tmp[1]]=nodeCount;
				tmp[1]=nodeCount;
				nodeCount++;
			}else{
				tmp[1]=nodeRelabel[tmp[1]];
			}
			fwrite(tmp,sizeof(u32),2,outputEdge);
		}
	}
//	printf("\n");
	fclose(outputEdge);
	fclose(outputNode);

	cudaFreeHost(nodeState);	
	cudaFree(edge_gpu);
	free(edge);
	free(nodeRelabel);
	cudaFree(idx_gpu);
	free(ances_stack);


		
	clock_t end=clock();
	double time_taken;
	time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("joinEdgeIndel: %.2fs elapsed for %u tips in GPU %d.\n",time_taken,ntip,runDevice);

    return 0;
}
