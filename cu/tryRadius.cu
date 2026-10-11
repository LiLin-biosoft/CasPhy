#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <ctime>

#define BLOCKSIZE 128
#define MAXMUTATE (1024*256)
#define MAXNEIGHBOR (1024*512)
#define NRESERVE 1024
#define WIDTH 32
#define LOG2WIDTH 5

#define CUDACALL(call) {                                    \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        printf("CUDA error at %s:%d code=%d(%s) \"%s\"\n",          \
               __FILE__, __LINE__, err, cudaGetErrorString(err), #call); \
        exit(1);                                                    \
    }                                                               \
}

typedef unsigned int u32;

__global__ void findNeighborKernel(u32* const mat,u32 const nrow,u32 const ntip,
							float radius,
							u32 const nMutate,u32* const mutateIdx,
							float* const weightRun,u32* nNeighbor){
	const u32 tid_local=threadIdx.x;
	const u32 thread_xidx=threadIdx.x+blockIdx.x*blockDim.x;
	float commonCount=0;
	if(thread_xidx<ntip){
		for(u32 i=0;i<nMutate;i++){
			u32 tmpIdx=mutateIdx[i];
			u32 xx=tmpIdx>>LOG2WIDTH;
			u32 yy=tmpIdx&(WIDTH-1);
			u32 tmp=mat[thread_xidx+xx*ntip];
			if(((tmp>>yy)&1)==1){
				commonCount+=weightRun[i];
			}
		}
	}
    u32 active = __activemask();
    u32 mask = __ballot_sync(active, commonCount >= radius);
    u32 nCount=__popc(mask);
    if((tid_local&31)==0){
   	atomicAdd(nNeighbor,nCount);
    }
}

int main(int argc, char *argv[]){
	setbuf(stdout, NULL); 
	clock_t start = clock();

	if(argc<3){
		printf("Error: no sufficient args.\n");
		return 0;
	} 

	float radius=atoi(argv[2]);
	u32 nReplicate = atoi(argv[3]);

	char input_file[256];
	strcpy(input_file,argv[1]);
	strcat(input_file,".bmat");
    u32 dim[2]; 
    FILE *file=fopen(input_file, "rb");
    if (file == NULL) {
    	printf("No bmat input.\n");
    	return 0;
  	} 
    fread(dim, sizeof(u32), 2, file);
    const u32 nrow=dim[0],ntip=dim[1];	
    const u32 nsite=nrow*WIDTH;
	u32* mat=(u32*)(malloc(sizeof(u32)*nrow*ntip));

	u32* mat_gpu;
	CUDACALL(cudaMalloc((void**)&mat_gpu,sizeof(u32)*nrow*ntip));
	for(u32 i=0;i<nrow;i++){
		fread(mat+i*ntip, sizeof(u32), ntip, file);
		CUDACALL(cudaMemcpy(mat_gpu+i*ntip,mat+i*ntip, sizeof(u32)*ntip, cudaMemcpyHostToDevice));
	}
	fclose(file);
	
	strcpy(input_file,argv[1]);
	strcat(input_file,".w");
    file=fopen(input_file, "rb");
    if (file == NULL) {
    	printf("No weight input.\n");
    	return 0;
  	} 
  	
  	double* weightTmp=(double*)(malloc(sizeof(double)*nsite));
	fread(weightTmp, sizeof(double), nsite, file);  
	fclose(file);	
	float* weight=(float*)(malloc(sizeof(float)*nsite));
	for(u32 i=0;i<nsite;i++){
		weight[i]=(float)(weightTmp[i]);
	}
	free(weightTmp);
	
	u32* mutateIdx;
	CUDACALL(cudaHostAlloc((void**)&mutateIdx,sizeof(u32)*MAXMUTATE,cudaHostAllocDefault));
	u32* mutateIdx_gpu;
	CUDACALL(cudaMalloc((void**)&mutateIdx_gpu,sizeof(u32)*MAXMUTATE));
	float* weightRun;
	CUDACALL(cudaHostAlloc((void**)&weightRun,sizeof(float)*MAXMUTATE,cudaHostAllocDefault));
	float* weightRun_gpu;
	CUDACALL(cudaMalloc((void**)&weightRun_gpu,sizeof(float)*MAXMUTATE));

	u32* idx=(u32*)(malloc(sizeof(u32)*nReplicate));
	u32* nNeighbor=(u32*)(malloc(sizeof(u32)*nReplicate));
	u32* nNeighbor_gpu;
	CUDACALL(cudaMalloc((void**)&nNeighbor_gpu,sizeof(u32)*nReplicate));
	CUDACALL(cudaMemset(nNeighbor_gpu,0,sizeof(u32)*nReplicate));

	char output_file[256];
	strcpy(output_file,argv[1]);
	
	std::srand((unsigned)std::time(nullptr));
	for(u32 i=0;i<nReplicate;i++){
		if((i&1023)==1023){
	      printf("=");
	    }
	    if((i&32767)==32767){
	      printf("\n");
	    }
		u32 rIdx=std::rand()%ntip;
		idx[i]=rIdx;
		
		u32 nMutate=0;
		for(u32 j=0;j<nrow;j++){
			u32 tmp=mat[rIdx+j*ntip];
			for(u32 k=0;k<WIDTH;k++){
				if(((tmp>>k)&1)==1){
					mutateIdx[nMutate]=j*WIDTH+k;
					weightRun[nMutate]=weight[j*WIDTH+k];
					nMutate++;
				}
			}
		}
		CUDACALL(cudaMemcpy(
			mutateIdx_gpu,mutateIdx,
			sizeof(u32)*nMutate,cudaMemcpyHostToDevice
		));
		CUDACALL(cudaMemcpy(
			weightRun_gpu,weightRun,
			sizeof(float)*nMutate,cudaMemcpyHostToDevice
		));
		
		findNeighborKernel<<<(ntip-1)/BLOCKSIZE+1,BLOCKSIZE>>>(
			mat_gpu,nrow,ntip,
			radius,
			nMutate,mutateIdx_gpu,
			weightRun_gpu,&(nNeighbor_gpu[i])
		);
	}
	printf("\n");
	CUDACALL(cudaMemcpy(
		nNeighbor,nNeighbor_gpu,sizeof(u32)*nReplicate,cudaMemcpyDeviceToHost
	));
	
	strcpy(output_file,argv[1]);
	strcat(output_file,".nNeighbor");
	FILE *output_nNeighbor=fopen(output_file,"wb");
	dim[0]=2,dim[1]=nReplicate;
	fwrite(dim,sizeof(u32),2,output_nNeighbor);
	fwrite(idx,sizeof(u32),nReplicate,output_nNeighbor);	
	fwrite(nNeighbor,sizeof(u32),nReplicate,output_nNeighbor);	
	fclose(output_nNeighbor);
	
	free(mat);
	free(idx);
	free(nNeighbor);
	CUDACALL(cudaFree(mat_gpu));
	CUDACALL(cudaFree(nNeighbor_gpu));
	CUDACALL(cudaFreeHost(mutateIdx));
	CUDACALL(cudaFree(mutateIdx_gpu));
	CUDACALL(cudaFreeHost(weightRun));
	CUDACALL(cudaFree(weightRun_gpu));
	
	
	clock_t end=clock();
	double time_taken;
	time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("radiusNeighbor: %.2fs elapsed for %u tips.\n",time_taken,ntip);
    return 0;
}
