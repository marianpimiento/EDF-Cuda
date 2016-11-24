//****librerias****

#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<time.h>

//****Variables globales****

int N=93, dimx=1920,dimy=2560,tam_imag=1920*2560;

//****Kernel: Función del device****

__global__ void Kernel(int *R_d, int *G_d, int *B_d, int *T_d, int *Rf, int *Gf, int *Bf, int d){
	int idx = threadIdx.x + blockIdx.x*blockDim.x;

	int tam_imag;
	tam_imag=1920*2560;
	if(idx<tam_imag)
		if(T_d[idx]==d){
			Rf[idx]=R_d[idx];
			Gf[idx]=G_d[idx];
			Bf[idx]=B_d[idx];
		}
}

//****Función main()****

int main(int argc,char* argv[]){

	//****declaración de variables para el host y device****

	int i, j, d, cont;
	int *R_h, *R_d, *G_h, *G_d, *B_h, *B_d, *T_d, *T_h;
	int *R, *G, *B;
	int *Rf, *Gf, *Bf;


	FILE *file, *Red, *Green, *Blue;
	FILE *ArchivoR, *ArchivoG, *ArchivoB;

	//****Leer archivo de la matriz topografica****

	file=fopen("Resultados/topo","r+");

	//******matriz R host y device********

	R_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&R_d, tam_imag*sizeof(int));

	//******matriz G host y device********

	G_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&G_d, tam_imag*sizeof(int));

	//******matriz B host y device********

	B_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&B_d, tam_imag*sizeof(int));

	//******matriz topografica device******

	T_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&T_d, tam_imag*sizeof(int));

	//**********matrices finales device*********

	cudaMalloc((void**)&Rf, tam_imag*sizeof(int));
	cudaMalloc((void**)&Gf, tam_imag*sizeof(int));
	cudaMalloc((void**)&Bf, tam_imag*sizeof(int));

	//******matrices Resultados finales******

	R=(int *)malloc(sizeof(int)*tam_imag);
	G=(int *)malloc(sizeof(int)*tam_imag);
	B=(int *)malloc(sizeof(int)*tam_imag);

	//***********cálculo del tiempo de procesamiento***********

	float t;
	clock_t tinicio, tfinal;
	tinicio=clock();

	//******matriz topografica***********
	cont=0;

	for(i=0;i<dimx;i++)
		for(j=0;j<dimy;j++){
			fscanf(file, "%d", &T_h[cont]);
			cont++;
		}
	fclose(file)

	//******operaciones*******

	for(d=1;d<=N;d++){

		//*******matriz R*********

		char ruta1[]="MiTesis/";
		sprintf(ruta1, "%s%d%s","RGB/",d,"/R");
		Red=fopen(ruta1,"r+");

		for(i=0;i<dimx*dimy;i++)
			fscanf(Red, "%d", &R_h[i]);
		fclose(Red);

		//*******matriz G*********

		char ruta2[]="MiTesis/";
		sprintf(ruta2, "%s%d%s","RGB/",d,"/G");
		Green=fopen(ruta2,"r+");

		for(i=0;i<dimx*dimy;i++)
			fscanf(Green, "%d", &G_h[i]);
		fclose(Green);

		//*******matriz B*********

		char ruta3[]="MiTesis/";
		sprintf(ruta3, "%s%d%s","RGB/",d,"/B");
		Blue=fopen(ruta3,"r+");

		for(i=0;i<dimx*dimy;i++)
			fscanf(Blue, "%d", &B_h[i]);
		fclose(Blue);

		//********copia de variables del Host al Device***********

		cudaMemcpy(R_d,R_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(G_d,G_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(B_d,B_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(T_d,T_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		//******llamado del kernel********

		Kernel<<<12288,400>>>(R_d, G_d, B_d, T_d, Rf, Gf, Bf, d);

	}//Fin for

	//copia de variables del Device al Host

	cudaMemcpy(R, Rf, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
	cudaMemcpy(G, Gf, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
	cudaMemcpy(B, Bf, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);

	//almacenamiento de las matrices resultantes*******

	ArchivoR=fopen("Resultados/R","w+");
	ArchivoG=fopen("Resultados/G","w+");
	ArchivoB=fopen("Resultados/B","w+");

	for(i=0;i<tam_imag;i++){
		if(i%dimy==0 && i!=0){
			fprintf(ArchivoR,"\n");
			fprintf(ArchivoG,"\n");
			fprintf(ArchivoB,"\n");
		}

		fprintf(ArchivoR,"%d ",R[i]);
		fprintf(ArchivoG,"%d ",G[i]);
		fprintf(ArchivoB,"%d ",B[i]);
	}

	fclose(ArchivoR);
	fclose(ArchivoG);
	fclose(ArchivoB);

	tfinal=clock();
	t = ((float)tfinal-(float)tinicio)/CLOCKS_PER_SEC;
	printf("tiempo de procesamiento:%6.3f s\n",t);

}//Fin función main()
