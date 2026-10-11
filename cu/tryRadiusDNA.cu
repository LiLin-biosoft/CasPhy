#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <ctime>

#define BLOCKSIZE 32
#define MAXNEIGHBOR (1024*256)
#define WIDTH 32
#define LOG2WIDTH 5
#define BASEPERINT 4
#define ROWCACHE 16384

typedef unsigned int u32;

__global__ void packMatKernel(u32* const matRaw,
                             const u32 nrow,const u32 nbinRaw,const u32 pitchRaw,
                             u32* const mat,const u32 nbin,const u32 pitch) {
	const u32 xlocal=threadIdx.x;
	const u32 rowId=blockIdx.x;
	const u32 stride=blockDim.x;

	u32* rowRaw_ptr=(u32*)((char*)(matRaw)+rowId*pitchRaw);
	u32* row_ptr=(u32*)((char*)(mat)+rowId*pitch);
	for(u32 i=xlocal; (i>>1)<nbin; i+=stride) {
		u32 value=0;
		if(i<nbinRaw) {
			u32 valueRaw=rowRaw_ptr[i];
			for (u32 j = 0; j < 4; ++j) {
				unsigned char nuc = (valueRaw >> (j * 8)) & 0x0F;
				value|=(nuc<<(j*4));
			}
		}
		u32 mask=__activemask();
		u32 value2=__shfl_down_sync(mask,value,1);
		u32 valuePacked=value|(value2<<16);
		if((xlocal&1)==0) {
			row_ptr[i>>1]=valuePacked;
		}
	}

}

__global__ void countNeighborKernel(u32* qry,u32* const mat,
									const u32 ntip,const u32 nbin,const u32 pitch,
                                   	u32* const nNeighbor,const u32 radius) {
	const u32 x_local=threadIdx.x;
	const u32 stride=blockDim.x;
	const u32 y_tid=threadIdx.y+blockIdx.y*blockDim.y;
	u32 commonCount=0;
	if(y_tid<ntip){
		for(u32 i=x_local; i<nbin; i+=stride) {
			u32* row_ptr=(u32*)((char*)(mat)+y_tid*pitch);
			u32 value=row_ptr[i]&qry[i];
			for (u32 j = 0; j < 8; ++j) {
				unsigned char nuc = (value >> (j * 4)) & 0x0F;
				if (nuc != 0x00 ) {
					commonCount ++;
				}
			}
		}
	}

    for (u32 offset = 16; offset > 0; offset>>=1) {
        commonCount += __shfl_down_sync(0xFFFFFFFF, commonCount, offset);
    }

	if(x_local==0&&y_tid<ntip){
		if(commonCount>=radius){
			atomicAdd(nNeighbor,1);
		}
	}
}

int main(int argc, char *argv[]) {
	setbuf(stdout, NULL);
	clock_t start = clock();

	if(argc < 4) {
		printf("Error: no sufficient args.\n");
		return 0;
	}

	u32 radius = atoi(argv[2]);
	u32 nReplicate = atoi(argv[3]);

	char input_file[256];
	strcpy(input_file, argv[1]);
	strcat(input_file, ".tbmat");

	u32 dim[2];
	FILE *file = fopen(input_file, "rb");
	if (file == NULL) {
		printf("No bmat input.\n");
		return 0;
	}
	fread(dim, sizeof(u32), 2, file);
	const u32 ntip = dim[0], nbinRaw = dim[1];

	u32* matRaw = (u32*)(malloc(sizeof(u32)*nbinRaw*ROWCACHE));
	u32* matRaw_gpu;
	size_t pitchRaw;
	cudaMallocPitch((void**)&matRaw_gpu, &pitchRaw, nbinRaw*sizeof(u32),ROWCACHE);
	u32* mat_gpu;
	u32 ncol=(nbinRaw-1)/2+1;
	size_t pitch;
	cudaMallocPitch((void**)&mat_gpu, &pitch, ncol*sizeof(u32),ntip);
	
	for(u32 i=0;i<ntip;i+=ROWCACHE){
		u32 nCache=ROWCACHE;
		if(ntip-i<nCache){
			nCache=ntip-i;
		}
		fread(matRaw, sizeof(u32), nbinRaw*nCache, file);
		cudaMemcpy2D(matRaw_gpu, pitchRaw, 
		             matRaw, nbinRaw * sizeof(u32),
		             nbinRaw * sizeof(u32), nCache,
		             cudaMemcpyHostToDevice);
		packMatKernel<<<nCache,32>>>(matRaw_gpu,nCache,nbinRaw,pitchRaw,
		                        (u32*)((char*)(mat_gpu)+i*pitch),ncol,pitch);
	}
	fclose(file);
	cudaFree(matRaw_gpu);
	free(matRaw);


	u32* nNeighbor_gpu;
	cudaMalloc((void**)&nNeighbor_gpu, sizeof(u32)*nReplicate);

	u32* idx=(u32*)(malloc(sizeof(u32)*nReplicate));
	u32* nNeighbor;
	cudaHostAlloc((void**)&nNeighbor, sizeof(u32)*nReplicate, cudaHostAllocDefault);
	cudaMemset(nNeighbor_gpu, 0, sizeof(u32)*nReplicate);

	const dim3 threads_square(BLOCKSIZE,BLOCKSIZE);
	const dim3 blocks_square(1,(ntip-1)/BLOCKSIZE+1);
	
	std::srand((unsigned)std::time(nullptr));
	for(u32 i = 0; i < nReplicate; i++) {
		u32 rIdx=std::rand()%ntip;
		idx[i]=rIdx;
		countNeighborKernel<<<blocks_square, threads_square>>>(
			(u32*)((char*)(mat_gpu)+rIdx*pitch),
		    mat_gpu, ntip, ncol,pitch,
		    nNeighbor_gpu + i,radius
		);
	}
	cudaMemcpy(nNeighbor, nNeighbor_gpu, sizeof(u32)*nReplicate, cudaMemcpyDeviceToHost);

	float avg=0;
	for(u32 i=0;i<nReplicate;++i){
		avg+=nNeighbor[i];
	}
	avg/=nReplicate;
	printf("Average nNeighbor: %.2f\n",avg);

	char output_file[256];
	strcpy(output_file, argv[1]);
	strcat(output_file, ".nNeighbor");
	FILE *output_nNeighbor = fopen(output_file, "wb");
	dim[0] = 2;
	dim[1] = nReplicate;
	fwrite(dim, sizeof(u32), 2, output_nNeighbor);
	fwrite(idx, sizeof(u32), nReplicate, output_nNeighbor);
	fwrite(nNeighbor, sizeof(u32), nReplicate, output_nNeighbor);
	fclose(output_nNeighbor);

	cudaFree(mat_gpu);
	cudaFree(nNeighbor_gpu);
	cudaFreeHost(nNeighbor);

	clock_t end = clock();
	double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("radiusNeighborDNA: %.2fs elapsed for %u tips.\n", time_taken, ntip);

	return 0;
}
