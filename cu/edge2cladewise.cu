#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>

#define BLOCKSIZE 128
#define NBLOCK 32

typedef unsigned int u32;

__global__ void findValueKernel(u32* const edge, u32 nedge,
                                u32 value, u32* idx) {
    __shared__ u32 threadFound[BLOCKSIZE];
    
    u32 localtid=threadIdx.x;
    u32 tid = threadIdx.x+blockIdx.x*blockDim.x;
    u32 stride = blockDim.x*gridDim.x;
    
    u32 lineIdx = nedge + 1;
    
    // each thread process several elements
    for(u32 i = tid; i < nedge; i += stride) {
        if(edge[i * 2] == value && i <= lineIdx) {
            lineIdx = i;
        }
    }
    
    // each thread save local minIdx to shared memory
    threadFound[localtid] = lineIdx;
    __syncthreads();
    
    // find local minIdx in each block
    for(u32 s = blockDim.x / 2; s > 0; s >>= 1) {
        if(localtid < s) {
            if(threadFound[localtid] > threadFound[localtid + s]) {
                threadFound[localtid] = threadFound[localtid + s];
            }
        }
        __syncthreads();
    }
    
    // save global minIdx to output
    if(localtid == 0) {
    	idx[blockIdx.x]=threadFound[0];
    }
}

__global__ void findMinKernel(u32* data,u32 size) {
	__shared__ u32 threadMin[BLOCKSIZE];
	
    u32 tid = threadIdx.x;
    u32 stride = blockDim.x;
    u32 tmp=1u<<32-1;
    for(u32 i = tid; i < size; i += stride) {
        if(data[i]<tmp) {
            tmp=data[i];
        }
    }
    threadMin[tid]=tmp;
    __syncthreads();
    
    for(u32 s = blockDim.x / 2; s > 0; s >>= 1) {
        if(tid < s) {
            if(threadMin[tid] > threadMin[tid + s]) {
                threadMin[tid] = threadMin[tid + s];
            }
        }
        __syncthreads();
    }
    if(tid == 0) {
    	data[0]=threadMin[0];
    }
}

int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);
    clock_t start = clock();
    
    FILE *file;
    char input_file[256];
    strcpy(input_file, argv[1]);
    strcat(input_file, ".edge"); //read mat,each row is one tip, node index is 0-based
    u32 dim[2]; 
    file = fopen(input_file, "rb");
    if (file == NULL) {
        printf("No edge input.\n");
        return 0;
    } 
    fread(dim, sizeof(u32), 2, file);
    const u32 nedge = dim[0];
    u32* edge;
    cudaHostAlloc((void**)&edge, sizeof(u32) * nedge * 2, cudaHostAllocDefault);
    u32* edge_gpu;
    cudaMalloc((void**)&edge_gpu, sizeof(u32) * nedge * 2);
    fread(edge, sizeof(u32), nedge * 2, file);
    cudaMemcpy(edge_gpu, edge, sizeof(u32) * nedge * 2, cudaMemcpyHostToDevice);
    fclose(file);

    const u32 ntip = atoi(argv[2]);    
    const u32 root = ntip;
    const u32 nnode = atoi(argv[3]);
    u32* nodeRelabel = (u32*)(malloc(sizeof(u32) * nnode));
    for(u32 i = 0; i < ntip; i++) {
        nodeRelabel[i] = i;
    }
    for(u32 i = ntip; i < nnode; i++) {
        nodeRelabel[i] = nnode;
    }
    u32 nodeCount = ntip;
    u32 idx;
    u32* idx_gpu;
    cudaMalloc((void**)&idx_gpu, sizeof(u32)*NBLOCK);    
    u32* ances_stack = (u32*)(malloc(sizeof(u32) * (nedge + 1)));
    ances_stack[0] = root;
    int stack_pointer = 0;

    FILE *outputEdge;
    char output_file[256];
    strcpy(output_file, argv[1]);
    strcat(output_file, ".edge2");
    outputEdge = fopen(output_file, "wb");
    dim[0] = nedge;
    dim[1] = 2;
    fwrite(dim, sizeof(u32), 2, outputEdge);
	u32 idxMax=1u<<32-1;
    while(stack_pointer > (-1)) {
        findValueKernel<<<NBLOCK, BLOCKSIZE>>>(edge_gpu, nedge, ances_stack[stack_pointer], idx_gpu); 
        findMinKernel<<<1,BLOCKSIZE>>>(idx_gpu,NBLOCK);
        cudaMemcpy(&idx, idx_gpu, sizeof(u32), cudaMemcpyDeviceToHost);
        if(idx >= nedge) {
            stack_pointer--;
        } else {
            cudaMemcpy(edge_gpu+idx*2,&idxMax,sizeof(u32),cudaMemcpyHostToDevice);
            stack_pointer++;
            ances_stack[stack_pointer] = edge[idx * 2 + 1];
            u32 tmp[2];
            tmp[0] = edge[idx * 2];
            tmp[1] = edge[idx * 2 + 1];
            
            if(nodeRelabel[tmp[0]] == nnode) {
                nodeRelabel[tmp[0]] = nodeCount;
                tmp[0] = nodeCount;
                nodeCount++;
            } else {
                tmp[0] = nodeRelabel[tmp[0]];
            }
            if(nodeRelabel[tmp[1]] == nnode) {
                nodeRelabel[tmp[1]] = nodeCount;
                tmp[1] = nodeCount;
                nodeCount++;
            } else {
                tmp[1] = nodeRelabel[tmp[1]];
            }
//            printf("write line %u,pointer: %d,root: %u\n",idx,stack_pointer,ances_stack[stack_pointer]);
            fwrite(tmp, sizeof(u32), 2, outputEdge);
        }
    }
    fclose(outputEdge);
    
    cudaFree(edge_gpu);
    cudaFreeHost(edge);
    free(nodeRelabel);
    cudaFree(idx_gpu);
    free(ances_stack);

    clock_t end=clock();
	double time_taken;
	time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("edge2cladewise: %.2fs elapsed for %u edges.\n",time_taken,nedge);   
    return 0; 
}
