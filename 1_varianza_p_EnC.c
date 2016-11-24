#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//************variables globales***************

int msk=3, dimx=36, dimy=36, tam_imag=36*36;

//*******************kernel********************



//*****************funcion main**********************

int main(int argc,char* argv[]){

	//***************declaracion de variables**************

	printf("\nInicia main");
	int i,j,m,cont,tam_B, init,fin;
	//init=atoi(argv[1]);
	//fin=atoi(argv[2]);
	init=1;
	fin=4;

	tam_B=(dimx+2)*(dimy+2);

	printf("\nCrea FILE");
	FILE *arch, *matrizG;

	int **A;
	int B[dimx+2][dimy+2];

	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	int *B_d, *B_h;
	float *var_d,*var_h;

	int d;
	printf("\nInicia for init-fin = %d - %d", init, fin);
	for(d=init;d<=fin;d++){

		//*******************declaracion de variables***************
		printf("\n\n***Entra for d=%d", d);
		
		B_h=(int *)malloc(sizeof(int)*tam_B);
		//cudaMalloc((void**)&B_d, tam_B*sizeof(int));
		var_h=(float *)malloc(sizeof(float)*tam_imag);
		//cudaMalloc((void**)&var_d,tam_imag*sizeof(float));
		
		A=(int **)malloc(sizeof(int)*dimx*dimy);
		for(i=0;i<dimx;i++)
			A[i]=(int*)malloc(sizeof(int)*dimy);


		//*****************calculo matriz B****************

		char ruta1[]="MiTesis/";
		sprintf(ruta1, "%s%d%s","RGB/",d,"/G"); 
		matrizG=fopen(ruta1,"r+");
		
		//int tem;
		for(i=0;i<dimx;i++){
			printf("\nEntra for con I = %d", i);
			for(j=0;j<dimy;j++){
				//printf("\nJ = %d", j);
				fscanf(matrizG, "%d", &A[i][j]);
				//printf("  A %d  ", A[i][j]);
				//fscanf(matrizG, "%d", &tem);
				//printf("  A %d  ", tem);
			}
		}
		fclose(matrizG);

		cont=0;
		for(i=0;i<dimx+2;i++){
			for(j=0;j<dimy+2;j++){
				B[i][j]=((i==0 || j==0 || i==dimx+1 || j==dimy+1) ? 0:A[i-1][j-1]);
				B_h[cont]=B[i][j];
				cont++;
			}

		}
/*
		//******************llamado de kernel*******************

		dim3 Grid(1);
		dim3 Block(36,36);

		cudaMemcpy(B_d,B_h,sizeof(int)*tam_B,cudaMemcpyHostToDevice);

		kernel<<<Grid,Block>>>(B_d,var_d);

		cudaMemcpy(var_h,var_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);


		//****************almacenamiento matriz de varianza**************

		char rutaV[]="VARIANZAS/";
		sprintf(rutaV, "%s%d", rutaV,d);
		arch=fopen(rutaV,"w+");

		for(m=0;m<tam_imag;m++){
			if(m%dimy==0 && m!=0){
				fprintf(arch,"\n");
			}
			fprintf(arch,"%f",var_h[m]);
		}
		fclose(arch);

		free(B_h);
		free(var_h);
		free(A);
		cudaFree(var_d);
		cudaFree(B_d);
		*/
		printf("\n***Fin for d=%d", d);
	}
	printf("\nFinaliza for init-fin = %d - %d", init, fin);


	// ???????????
/*
	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento para calcular varianzas de %d matrices: %6.3fs\n",fin-init+1,t);
*/
	return 0;

}//FIN funcion main()

