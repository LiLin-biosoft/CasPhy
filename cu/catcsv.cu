#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
	
	if(argc<3){
		printf("Error: no sufficient args.\n");
		return 0;
	} 
	
    FILE *file;
    file = fopen(argv[1], "rb");
    if(file==NULL){
    	printf("Error: input file.\n");
    	return 0;
	}
    unsigned int dim[2]; 
    fread(dim, sizeof(dim), 1, file);
    const unsigned int nrow=dim[0],ncol=dim[1];
    const unsigned int offset_i=atoi(argv[2]);
    const unsigned int offset_j=atoi(argv[3]);
    printf("read matrix with %d rows and %d cols\n",nrow,ncol);
    printf("print from %dth rows and %dth cols\n",offset_i,offset_j);

    unsigned int* mat=(unsigned int*)(malloc(ncol*sizeof(unsigned int)));
    for(unsigned int i=0;i<offset_i;i++){
    	fseek(file, sizeof(unsigned int)*ncol, SEEK_CUR);
	}

	for(unsigned int i=0;i<16&&i+offset_i<nrow;i++){
		fread(mat,sizeof(unsigned int),ncol,file);
    	for(unsigned int j=0;j<16&&j+offset_j<ncol;j++){
    		printf("%d,",mat[j+offset_j]);
		}
		printf("\n");
	}

	fclose(file);
    return 0;
}
