#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>

#define BLOCKSIZE 32
#define CROSSSIZE 8
#define CACHESIZE (1024*2)

#define CUDACALL(call) {                                    \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        printf("CUDA error at %s:%d code=%d(%s) \"%s\"\n",          \
               __FILE__, __LINE__, err, cudaGetErrorString(err), #call); \
        exit(1);                                                    \
    }                                                               \
}

//TransferDisMulti

typedef unsigned int u32;

__global__ void crossXorKernel(u32* const mat_x,u32* const mat_y,
							  u32* const xmin,u32* const ymin,
							  u32 const nrow_x,u32 const nrow_y,u32 const ncol,u32 const ntip){
	__shared__ u32 threadSum1[CROSSSIZE*CROSSSIZE];
	__shared__ u32 threadSum2[CROSSSIZE*CROSSSIZE];
	__shared__ u32 threadSumComp1[CROSSSIZE*CROSSSIZE];
	__shared__ u32 threadSumComp2[CROSSSIZE*CROSSSIZE];
	
	__shared__ u32 xCache[CROSSSIZE*CROSSSIZE];
	__shared__ u32 yCache[CROSSSIZE*CROSSSIZE];
	
	const u32 xlocal=threadIdx.x;
	const u32 ylocal=threadIdx.y;
	const u32 xtid=threadIdx.x+blockIdx.x*blockDim.x;
	const u32 ytid=threadIdx.y+blockIdx.y*blockDim.y;
	const u32 xtidread=threadIdx.y+blockIdx.x*blockDim.x;
	
	threadSum1[ylocal*CROSSSIZE+xlocal]=ncol*32;
	threadSum2[ylocal*CROSSSIZE+xlocal]=ncol*32;
	threadSumComp1[ylocal*CROSSSIZE+xlocal]=ncol*32;
	threadSumComp2[ylocal*CROSSSIZE+xlocal]=ncol*32;
	__syncthreads();
	
	u32 count=0;
	u32 countComp=0;
	for(u32 i=0;i<ncol;i+=CROSSSIZE){
		if(ytid<nrow_y&&i+xlocal<ncol){
			yCache[ylocal*CROSSSIZE+xlocal]=mat_y[ytid*ncol+i+xlocal];
		}
		if(xtidread<nrow_x&&i+xlocal<ncol){
			xCache[ylocal*CROSSSIZE+xlocal]=mat_x[xtidread*ncol+i+xlocal];
		}
		__syncthreads();

		for(u32 j=0;j<CROSSSIZE&&i+j<ncol;j++){
			if(xtid<nrow_x&&ytid<nrow_y){
				const u32 tmp=xCache[xlocal*CROSSSIZE+j]^yCache[ylocal*CROSSSIZE+j];
				const u32 tmpComp=xCache[xlocal*CROSSSIZE+j]^(~yCache[ylocal*CROSSSIZE+j]);
				count+=__popc(tmp);
				countComp+=__popc(tmpComp);
			}			
		}
		__syncthreads();
	}
	if(xtid<nrow_x&&ytid<nrow_y){
		threadSum1[ylocal*CROSSSIZE+xlocal]=count;
		threadSum2[ylocal+xlocal*CROSSSIZE]=count;
		threadSumComp1[ylocal*CROSSSIZE+xlocal]=countComp;
		threadSumComp2[ylocal+xlocal*CROSSSIZE]=countComp;
	}
	__syncthreads();
	
    for(u32 s=blockDim.y/2;s>0;s>>=1){
        if(ylocal<s){
            if(threadSum1[ylocal*CROSSSIZE+xlocal]>threadSum1[(ylocal+s)*CROSSSIZE+xlocal]){
				threadSum1[ylocal*CROSSSIZE+xlocal]=threadSum1[(ylocal+s)*CROSSSIZE+xlocal];
            }
            if(threadSumComp1[ylocal*CROSSSIZE+xlocal]>threadSumComp1[(ylocal+s)*CROSSSIZE+xlocal]){
				threadSumComp1[ylocal*CROSSSIZE+xlocal]=threadSumComp1[(ylocal+s)*CROSSSIZE+xlocal];
            }
            if(threadSum2[ylocal*CROSSSIZE+xlocal]>threadSum2[(ylocal+s)*CROSSSIZE+xlocal]){
				threadSum2[ylocal*CROSSSIZE+xlocal]=threadSum2[(ylocal+s)*CROSSSIZE+xlocal];
            }
            if(threadSumComp2[ylocal*CROSSSIZE+xlocal]>threadSumComp2[(ylocal+s)*CROSSSIZE+xlocal]){
				threadSumComp2[ylocal*CROSSSIZE+xlocal]=threadSumComp2[(ylocal+s)*CROSSSIZE+xlocal];
            }
        }
        __syncthreads();
    }
    if(ylocal==0){
    	u32 tmp1=threadSum1[xlocal];
    	u32 tmp2=threadSumComp1[xlocal]-(ncol*32-ntip);
    	if(tmp1>tmp2){
    		tmp1=tmp2;
		}
	    atomicMin(&(xmin[xtid]),tmp1);  
		tmp1=threadSum2[xlocal];
		tmp2=threadSumComp2[xlocal]-(ncol*32-ntip);
		if(tmp1>tmp2){
			tmp1=tmp2;		
		}
		atomicMin(&(ymin[threadIdx.x+blockIdx.y*blockDim.y]),tmp1);
	}
}

__global__ void rowSizeKernel(u32* const mat,u32* const size,
							 u32 const nrow,u32 const ncol,u32 const ntip){
	const u32 thread_xidx=threadIdx.x+blockIdx.x*blockDim.x;
	if(thread_xidx<nrow){
		u32 count=0;
		for(u32 i=0;i<ncol;i++){
			const u32 tmp=mat[thread_xidx*ncol+i];
			count+=__popc(tmp);
		}
		if(count>ntip-count){
			count=ntip-count;
		}
		size[thread_xidx]=count;
	}
}

__global__ void fillvalue(u32* const data,const u32 value,const u32 size){
	const u32 idx=(blockIdx.x*blockDim.x)+threadIdx.x;
	if(idx<size){
		data[idx]=value;
	}
}

int main(int argc, char *argv[]){
	setbuf(stdout, NULL); 
	if(argc<3){
		printf("Error: no sufficient args.\n");
		return 0;
	} 
	const int gpuCount=2;
	
	u32 ntip=atoi(argv[4]);
	FILE *file_y;
    u32 dim[2]; 
    file_y = fopen(argv[1], "rb");
    fread(dim, sizeof(u32), 2, file_y);
    const u32 nrow_y=dim[0],ncol=dim[1];
    FILE *file_x;
    file_x = fopen(argv[2], "rb");
    fread(dim, sizeof(u32), 2, file_x);
    const u32 nrow_x=dim[0];   

    cudaStream_t* s_run=(cudaStream_t*)(malloc(sizeof(cudaStream_t)*gpuCount));
	for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		CUDACALL(cudaStreamCreate(&s_run[gpuIdx]));
	}

	CUDACALL(cudaSetDevice(0));
	u32* ysize_gpu;
	CUDACALL(cudaMalloc((void**)&ysize_gpu,sizeof(u32)*nrow_y));
	CUDACALL(cudaMemsetAsync(ysize_gpu,0,sizeof(u32)*nrow_y,s_run[0]));
	
	u32** buf=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32** buf_gpu=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32** mat_y_gpu=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32** xsize_gpu=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32** ymin_gpu=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32** xmin_gpu=(u32**)(malloc(sizeof(u32*)*gpuCount));
	u32* xsize;
	CUDACALL(cudaHostAlloc((void**)&xsize,sizeof(u32)*nrow_x,cudaHostAllocDefault));
	
	for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		CUDACALL(cudaHostAlloc(
			(void**)&(buf[gpuIdx]),
			sizeof(u32)*CACHESIZE*ncol,
			cudaHostAllocDefault
		));
		CUDACALL(cudaMalloc((void**)&(buf_gpu[gpuIdx]), sizeof(u32)*CACHESIZE*ncol));
		CUDACALL(cudaMalloc((void**)&(mat_y_gpu[gpuIdx]), sizeof(u32)*CACHESIZE*ncol));
		CUDACALL(cudaMalloc((void**)&xsize_gpu[gpuIdx],sizeof(u32)*CACHESIZE));
		CUDACALL(cudaMalloc((void**)&ymin_gpu[gpuIdx],sizeof(u32)*nrow_y));
		CUDACALL(cudaMalloc((void**)&xmin_gpu[gpuIdx],sizeof(u32)*nrow_x));
		CUDACALL(cudaMemsetAsync(xsize_gpu[gpuIdx],0,sizeof(u32)*CACHESIZE,s_run[gpuIdx]));
		fillvalue<<<(nrow_y-1)/BLOCKSIZE+1,BLOCKSIZE,0,s_run[gpuIdx]>>>(ymin_gpu[gpuIdx],ncol*32,nrow_y);
		fillvalue<<<(nrow_x-1)/BLOCKSIZE+1,BLOCKSIZE,0,s_run[gpuIdx]>>>(xmin_gpu[gpuIdx],ncol*32,nrow_x);
	}
	
	for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		cudaStreamSynchronize(s_run[gpuIdx]);
	}

	printf("calculate crossXOR between 1st mat with %d rows and 2nd mat with %d rows.\n",nrow_y,nrow_x);
	clock_t start = clock();
	for(int i=0;i<nrow_y;i+=CACHESIZE){
		u32 length_y=CACHESIZE;
		if(length_y>nrow_y-i){
			length_y=nrow_y-i;
		}
		fread(buf[0], sizeof(u32), length_y*ncol, file_y);
		for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
			CUDACALL(cudaSetDevice(gpuIdx));
			CUDACALL(cudaMemcpyAsync(
				mat_y_gpu[gpuIdx],buf[0], sizeof(u32)*length_y*ncol,
				cudaMemcpyHostToDevice,s_run[gpuIdx]
			));			
		}
		CUDACALL(cudaSetDevice(0));
		rowSizeKernel<<<(length_y-1)/BLOCKSIZE+1,BLOCKSIZE,0,s_run[0]>>>(
			mat_y_gpu[0],ysize_gpu+i,length_y,ncol,ntip
		);
		fseek(file_x, sizeof(u32)*2, SEEK_SET);
		u32 gpuIdx=0;
		
		for(int j=0;j<nrow_x;j+=CACHESIZE){
			u32 length_x=CACHESIZE;
			if(length_x>nrow_x-j){
				length_x=nrow_x-j;
			}
			CUDACALL(cudaSetDevice(gpuIdx));
			CUDACALL(cudaStreamSynchronize(s_run[gpuIdx]));
			fread(buf[gpuIdx], sizeof(u32), length_x*ncol, file_x);
			CUDACALL(cudaMemcpyAsync(
				buf_gpu[gpuIdx],buf[gpuIdx], sizeof(u32)*length_x*ncol,
				cudaMemcpyHostToDevice,s_run[gpuIdx]
			));
			printf("=");
			if(i==0){
				rowSizeKernel<<<(length_x-1)/BLOCKSIZE+1,BLOCKSIZE,0,s_run[gpuIdx]>>>(
					buf_gpu[gpuIdx],xsize_gpu[gpuIdx],length_x,ncol,ntip
				);
				CUDACALL(cudaMemcpyAsync(
					xsize+j,
					xsize_gpu[gpuIdx],
					sizeof(u32)*length_x,
					cudaMemcpyDeviceToHost,
					s_run[gpuIdx]
				));
			}

			const dim3 blocks_square((length_x-1)/CROSSSIZE+1,(length_y-1)/CROSSSIZE+1);
			const dim3 threads_square(CROSSSIZE,CROSSSIZE);	
			crossXorKernel<<<blocks_square,threads_square,0,s_run[gpuIdx]>>>(
				buf_gpu[gpuIdx],mat_y_gpu[gpuIdx],
				xmin_gpu[gpuIdx]+j,ymin_gpu[gpuIdx]+i,
				length_x,length_y,
				ncol,ntip
			);
			gpuIdx=(gpuIdx+1)%gpuCount;
		}
		
		for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
			CUDACALL(cudaSetDevice(gpuIdx));
			cudaStreamSynchronize(s_run[gpuIdx]);
		}
		printf("\n");
	}

	u32* xmin;
	CUDACALL(cudaHostAlloc((void**)&xmin,sizeof(u32)*nrow_x,cudaHostAllocDefault));
	u32 gpuIdx=0;
	for(u32 j=0;j<nrow_x;j+=CACHESIZE){
		u32 length_x=CACHESIZE;
		if(length_x>nrow_x-j){
			length_x=nrow_x-j;
		}
		CUDACALL(cudaSetDevice(gpuIdx));
		CUDACALL(cudaMemcpyAsync(
			xmin+j,
			xmin_gpu[gpuIdx]+j,
			sizeof(u32)*length_x,
			cudaMemcpyDeviceToHost,
			s_run[gpuIdx]
		));		
		gpuIdx=(gpuIdx+1)%gpuCount;
	}

	for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		cudaStreamSynchronize(s_run[gpuIdx]);
	}

	CUDACALL(cudaSetDevice(0));
	u32* ysize=(u32*)(malloc(sizeof(u32)*nrow_y));
	CUDACALL(cudaMemcpy(ysize, ysize_gpu, sizeof(u32)*nrow_y, cudaMemcpyDeviceToHost));
	CUDACALL(cudaFree(ysize_gpu)); 
	
	u32* ymin=(u32*)(malloc(sizeof(u32)*nrow_y));
	CUDACALL(cudaMemcpy(ymin, ymin_gpu[0], sizeof(u32)*nrow_y, cudaMemcpyDeviceToHost));
	u32* yminLocal=(u32*)(malloc(sizeof(u32)*nrow_y));
	for(u32 gpuIdx=1;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		CUDACALL(cudaMemcpy(yminLocal, ymin_gpu[gpuIdx], sizeof(u32)*nrow_y, cudaMemcpyDeviceToHost));
		for(u32 i=0;i<nrow_y;++i){
			if(yminLocal[i]<ymin[i]){
				ymin[i]=yminLocal[i];
			}
		}
	}
	free(yminLocal);
	
	for(u32 gpuIdx=0;gpuIdx<gpuCount;++gpuIdx){
		CUDACALL(cudaSetDevice(gpuIdx));
		CUDACALL(cudaStreamDestroy(s_run[gpuIdx]));
		CUDACALL(cudaFreeHost(buf[gpuIdx]));
		CUDACALL(cudaFree(buf_gpu[gpuIdx]));
		CUDACALL(cudaFree(ymin_gpu[gpuIdx])); 
		CUDACALL(cudaFree(xmin_gpu[gpuIdx])); 
		CUDACALL(cudaFree(xsize_gpu[gpuIdx])); 
	}
	
	clock_t end=clock();
	double time_taken;
	time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Running time on GPU: %.5fs\n",time_taken);
	
	
	fclose(file_y);
	fclose(file_x);
	
	FILE *file;
	char output_file[256];
	strcpy(output_file,argv[3]);
	strcat(output_file,"_ymin.bin");
	file=fopen(output_file,"wb");
	dim[0]=2,dim[1]=nrow_y;
	fwrite(dim,sizeof(u32),2,file);
	fwrite(ymin,sizeof(u32),nrow_y,file);
	fwrite(ysize,sizeof(u32),nrow_y,file);
	fclose(file);
	free(ymin);
	free(ysize); 
	
	strcpy(output_file,argv[3]);
	strcat(output_file,"_xmin.bin");
	file=fopen(output_file,"wb");
	dim[0]=2,dim[1]=nrow_x;
	fwrite(dim,sizeof(u32),2,file);
	fwrite(xmin,sizeof(u32),nrow_x,file);
	fwrite(xsize,sizeof(u32),nrow_x,file);
	fclose(file);	
	CUDACALL(cudaFreeHost(xmin));
	CUDACALL(cudaFreeHost(xsize));

    return 0;
}
