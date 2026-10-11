#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#define BLOCKSIZE 32
#define ROWCACHE 16384

typedef unsigned int u32;

__global__ void multipleKernel(u32* data, u32 size){
	u32 tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < size) {
        data[tid] *= 2;
    }
}

/* __global__ void sampleMatKernel(u32* const data,const u32 nrow,const u32 nbin,const u32 pitch,
								u32* const output, const u32 output_pitch,
								u32* const idx, const u32 nIdx){
	const u32 tid=threadIdx.x;
	const u32 bid=blockIdx.x;
	u32* data_row=(u32*)((char*)(data)+bid*pitch);
	u32* output_row=(u32*)((char*)(output)+bid*output_pitch);
	for(u32 i=tid;i<nbin;i+=BLOCKSIZE){
		output_row[i]=data_row[i];
	}
} */

__global__ void sampleMatKernel(u32* const __restrict__ data,const u32 nbin,const u32 pitch,
								u32* const __restrict__ output, const u32 output_pitch,
								u32* const __restrict__ idx, const u32 nIdx){
	const u32 tid=threadIdx.x;
	const u32 bid=blockIdx.x;
	u32* data_row=(u32*)((char*)(data)+bid*pitch);
	u32* output_row=(u32*)((char*)(output)+bid*output_pitch);
	for(u32 i=tid;i*4<nIdx;i+=BLOCKSIZE){
		u32 value=0;
		for(u32 j=0;(j<4)&&(i*4+j<nIdx);++j){
			u32 tmp=idx[i*4+j];
			u32 xx=tmp>>2;
			u32 yy=tmp&3;
			if(xx<nbin){
				u32 tmpInput=(data_row[xx]>>(yy*8))&0xff;
				value|=tmpInput<<(j*8);	
			}
		}
		output_row[i]=value;
	}
}

__host__ thrust::device_vector<u32> generateIdx(u32 size,u32 min_val, u32 max_val,u32 seed){
	std::mt19937 gen(seed);
	std::uniform_int_distribution<u32> dist(min_val, max_val);
	thrust::host_vector<u32> h_data(size);
    for (u32 i = 0; i < size; ++i) {
        h_data[i] = dist(gen);
    }
	thrust::device_vector<u32> d_data = h_data;
	thrust::sort(d_data.begin(), d_data.end());
	return d_data;
}

int main(int argc, char *argv[]) {
	setbuf(stdout, NULL);
	clock_t start = clock();

	if(argc < 3) {
		printf("Error: no sufficient args.\n");
		return 0;
	}

	u32 seed=atoi(argv[2]);
	char bsIdx[20];
	sprintf(bsIdx,"B%d",seed);

	char input_file[256];
	strcpy(input_file, argv[1]);
	strcat(input_file, ".bin");
	char output_file[256];
	strcpy(output_file, argv[1]);
	strcat(output_file, bsIdx);
	strcat(output_file, ".bin");
	char idx_file[256];
	strcpy(idx_file, argv[1]);
	strcat(idx_file, bsIdx);
	strcat(idx_file, ".idx");

	u32 dim[2];
	FILE *input = fopen(input_file, "rb");
	if (input == NULL) {
		printf("No input.\n");
		return 0;
	}
	fread(dim, sizeof(u32), 2, input);
	const u32 ntip = dim[0], nbin = dim[1];
	FILE *output = fopen(output_file, "wb");
	fwrite(dim, sizeof(u32), 2, output);
	
	u32* cache = (u32*)(malloc(sizeof(u32)*nbin*ROWCACHE));
	u32* cache_gpu;
	size_t pitch;
	cudaMallocPitch((void**)&cache_gpu, &pitch, nbin*sizeof(u32),ROWCACHE);
	u32* output_gpu;
	size_t output_pitch;
	cudaMallocPitch((void**)&output_gpu, &output_pitch, nbin*sizeof(u32),ROWCACHE);
	
	thrust::device_vector<u32> idx_gpu=generateIdx(nbin*4, 0, nbin*4-1,seed);
	dim[0]=1;dim[1]=nbin*4;
	thrust::host_vector<u32> idx = idx_gpu;
	FILE *idx_output = fopen(idx_file, "wb");
	fwrite(dim, sizeof(u32), 2, idx_output);
	fwrite(idx.data(),sizeof(u32),nbin*4,idx_output);
	fclose(idx_output);
	
	for(u32 i=0;i<ntip;i+=ROWCACHE){
		u32 nCache=ROWCACHE;
		if(ntip-i<nCache){
			nCache=ntip-i;
		}
		fread(cache, sizeof(u32), nbin*nCache, input);
		cudaMemcpy2D(cache_gpu, pitch, 
		             cache, nbin * sizeof(u32),
		             nbin * sizeof(u32), nCache,
		             cudaMemcpyHostToDevice);
		sampleMatKernel<<<nCache,32>>>(cache_gpu,nbin,pitch,
									   output_gpu,output_pitch,
									   thrust::raw_pointer_cast(idx_gpu.data()),nbin*4);
		cudaMemcpy2D(cache, nbin * sizeof(u32), 
		             output_gpu, output_pitch,
		             nbin * sizeof(u32), nCache,
		             cudaMemcpyDeviceToHost);
		fwrite(cache, sizeof(u32), nbin*nCache, output);
	}
	fclose(input);
	fclose(output);
	
	free(cache);
	cudaFree(cache_gpu);
	cudaFree(output_gpu);

	clock_t end = clock();
	double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("Felsenstein sampling: %.2fs elapsed for %u tips.\n", time_taken, ntip);

	return 0;
}
