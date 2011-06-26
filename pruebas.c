#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>


void main()
{
    //2D array
    double **array;
    
    double *X;
    int N, M, i, j;
   
    printf("EXAMPLE OF DYNAMIC MEMORY ALLOCATION\n");
    printf("Array lenght (columns): ");
    scanf("%d",&N);
    printf("Array lenght (raws): ");
    scanf("%d",&M);
    getchar();

    array = (double **)malloc (M * sizeof(double *));
    for ( i=0; i<M; i++ )
	array[i] = (double *)malloc ( N * sizeof(double));


    for( i=0; i<M; i++){
	for( j=0; j<N; j++){
	    array[i][j] = i*j*3.14159265358979 + i + j;
	    printf("X[%d][%d]: %lf\t",i,j,array[i][j]);}
	printf("\n");}
	
    for ( i=0; i<M; i++ )
      free (array[i]);
    free(array);

}