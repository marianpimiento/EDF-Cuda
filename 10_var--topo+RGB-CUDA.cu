
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h>

//************variables globales***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;

// [i][j] = i*dimy+j

//************** Kernel CUDA *********************
__global__ void Varianza (int *G_d, float *var_d){

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int id = idx + idy*blockDim.x*gridDim.x;

	int M_d[9], i, dimx=1040, dimy=1388, tam_imag=1388*1040, msk=3;
	float X=0.f,Xprom=0.f,Y=0.f;
	var_d[id]=0;
	//printf("prueba\n");

	if(id<tam_imag){
		//M_d[0]=((i<1 || j<1) ? 0:A[i-1][j-1]);
		/*
		M_d[0]=((idx<1 || idy<1) ? 0:G_d[(idx-1)+(idy-1)*blockDim.x*gridDim.x]);
		M_d[1]=((idx<1) ? 0:G_d[(idx-1)+(idy)*blockDim.x*gridDim.x]);
		M_d[2]=((idx<1 || idy>dimy-2) ? 0:G_d[(idx-1)+(idy+1)*blockDim.x*gridDim.x]);
		M_d[3]=((idy<1) ? 0:G_d[(idx)+(idy-1)*blockDim.x*gridDim.x]);
		M_d[4]=G_d[(idx)+(idy)*blockDim.x*gridDim.x];
		M_d[5]=((idy>dimy-2) ? 0:G_d[(idx)+(idy+1)*blockDim.x*gridDim.x]);
		M_d[6]=((idx>dimx-2 || idy<1) ? 0:G_d[(idx+1)+(idy-1)*blockDim.x*gridDim.x]);
		M_d[7]=((idx>dimx-2) ? 0:G_d[(idx+1)+(idy)*blockDim.x*gridDim.x]);
		M_d[8]=((idx>dimx-2 || idy>dimy-1) ? 0:G_d[(idx+1)+(idy+1)*blockDim.x*gridDim.x]);
 		*/
		
		if (idx==0 || idy==0){
			M_d[0]=0;
		}else{
			M_d[0]=G_d[id-1-dimy];
		}
/*
		if ((idx==0)){
			M_d[1]=0;
		}else{
			M_d[1]=G_d[id-dimy];
			//M_d[1]=8;
		}
/*
		if (idx==0 || idy==dimy){
			M_d[2]=0;
		}else{
			M_d[2]=G_d[id+1-dimy];
		}
*/
		if (idy==0){
			M_d[3]=0;
		}else{
			M_d[3]=G_d[id-1];
		}

		M_d[4]=G_d[id];

		if (idy==dimy){
			M_d[5]=0;
		}else{
			M_d[5]=G_d[id+1];
		}
/*
		if (id==dimx || idy==0){
			M_d[6]=0;
		}else{
			M_d[6]=G_d[id-1+dimy];
		}
*//*
		if (idx==dimx){
			M_d[7]=0;
		}else{
			M_d[7]=G_d[id+dimy];
		}
*//*
		if (idx==dimx || idy==dimy){
			M_d[8]=0;
		}else{
			M_d[8]=G_d[id+1+dimy];
		}
*/
		
		//M_d[0]=1;
		M_d[1]=5;
		M_d[2]=8;
		//M_d[3]=1;
		//M_d[4]=1;
		//M_d[5]=1;
		M_d[6]=2;
		M_d[7]=5;
		M_d[8]=4;
		
		for(i=0;i<msk*msk;i++)
			X+=M_d[i];

		Xprom=((float)X)/(msk*msk);

		for(i=0;i<msk*msk;i++)
			Y+=(Xprom-M_d[i])*(Xprom-M_d[i]);
		
		var_d[id]=Y/(msk*msk);

	}
}


//*****************Funcion main**********************
int main(int argc,char* argv[]){

	//***************Declaracion de variables**************
	int i,j,init,fin,d;
	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	//init=1;
	//fin=328;

	FILE *matrizR, *matrizG, *matrizB;
	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	int *topof_h, *R_h, *G_h, *B_h, *Rf_h, *Gf_h, *Bf_h;
	float *max_h, *var_h;

	int *topof_d, *R_d, *G_d, *B_d, *Rf_d, *Gf_d, *Bf_d;
	float *max_d, *var_d;
	
	//************Inicializacion de variables en el host y en el device *************** 
	
	/* // Declaracion tipo MATRIZ
	max_h=(float **)malloc(sizeof(float)*dimx);
	topof_h=(int **)malloc(sizeof(int)*dimx);
	R_h=(int **)malloc(sizeof(int)*dimx);
	G_h=(int **)malloc(sizeof(int)*dimx);
	B_h=(int **)malloc(sizeof(int)*dimx);
	Rf_h=(int **)malloc(sizeof(int)*dimx);
	Gf_h=(int **)malloc(sizeof(int)*dimx);
	Bf_h=(int **)malloc(sizeof(int)*dimx);
	for(i=0;i<dimx;i++){
		max_h[i]=(float*)malloc(sizeof(float)*dimy);
		topof_h[i]=(int*)malloc(sizeof(int)*dimy);
		R_h[i]=(int*)malloc(sizeof(int)*dimy);
		G_h[i]=(int*)malloc(sizeof(int)*dimy);
		B_h[i]=(int*)malloc(sizeof(int)*dimy);
		Rf_h[i]=(int*)malloc(sizeof(int)*dimy);
		Gf_h[i]=(int*)malloc(sizeof(int)*dimy);
		Bf_h[i]=(int*)malloc(sizeof(int)*dimy);
	}
	var_h=(float *)malloc(sizeof(float)*tam_imag);
	*/

	R_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&R_d, tam_imag*sizeof(int));
	G_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&G_d, tam_imag*sizeof(int));
	B_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&B_d, tam_imag*sizeof(int));
	Rf_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&Rf_d, tam_imag*sizeof(int));
	Gf_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&Gf_d, tam_imag*sizeof(int));
	Bf_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&Bf_d, tam_imag*sizeof(int));
	topof_h=(int *)malloc(sizeof(int)*tam_imag);
	cudaMalloc((void**)&topof_d, tam_imag*sizeof(int));

	max_h=(float *)malloc(sizeof(float)*tam_imag);
	cudaMalloc((void**)&max_d, tam_imag*sizeof(float));
	var_h=(float *)malloc(sizeof(float)*tam_imag);
	cudaMalloc((void**)&var_d,tam_imag*sizeof(float));

	
	//*************** For cálculo EDF ****************
	for(d=init;d<=fin;d++){

		printf("d=%d \n", d);
		
		//*****************Lecura de matrices RGB en el host****************

		char rutaR[]="";
		sprintf(rutaR, "%s%d%s","RGB/",d,"/R"); 
		matrizR=fopen(rutaR,"r+");

		char rutaG[]="";
		sprintf(rutaG, "%s%d%s","RGB/",d,"/G"); 
		matrizG=fopen(rutaG,"r+");

		char rutaB[]="";
		sprintf(rutaB, "%s%d%s","RGB/",d,"/B"); 
		matrizB=fopen(rutaB,"r+");

		for(i=0;i<dimx;i++){
			for(j=0;j<dimy;j++){
				fscanf(matrizR, "%d", &R_h[i*dimy+j]);
				fscanf(matrizG, "%d", &G_h[i*dimy+j]); 
				fscanf(matrizB, "%d", &B_h[i*dimy+j]); 
			}
		}
		fclose(matrizR);
		fclose(matrizG);
		fclose(matrizB);

		//***************** Kernel Varianza *******************

		cudaMemcpy(G_d,G_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		dim3 Grid(347,20);
		dim3 Block(13,16);

		Varianza<<<Grid,Block>>>(B_d,var_d);

		printf("Despues de kernel \n");

		cudaMemcpy(var_h,var_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);
		
		printf("Despues de resultado a host \n");
		//***************** Kernel Varianza *******************
		/*
		cudaMemcpy(R_d,R_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(G_d,G_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(B_d,B_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		cudaMemcpy(Rf_d,Rf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Gf_d,Gf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Bf_d,Bf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		cudaMemcpy(topof_d,topof_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(max_d,max_h,sizeof(float)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(var_d,var_h,sizeof(float)*tam_imag,cudaMemcpyHostToDevice);

		dim3 Grid(347,20);
		dim3 Block(13,16);

		TopoRGB<<<Grid,Block>>>(R_d,G_d,B_d,Rf_d,Gf_d,Bf_d,topof_d,max_d,var_d);

		cudaMemcpy(Rf_h,Rf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Gf_h,Gf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Bf_h,Bf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);

		cudaMemcpy(topof_h,topof_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(max_h,max_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);

		*/
		//*********************Calculo de TODO ********************
		
	} //Finaliza For cálculo EDF
	printf("***Sale del for \n");

	/*
	// ***************** Generacion de archivos de resultados ************************
	FILE *archTopo, *archR, *archG, *archB;
	archTopo=fopen("Resultados/topos10","w+");
	archR=fopen("Resultados/R10","w+");
	archG=fopen("Resultados/G10","w+");
	archB=fopen("Resultados/B10","w+");
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(archTopo,"%d ",topof_h[i*dimy+j]);
			fprintf(archR,"%d ",Rf_h[i*dimy+j]);
			fprintf(archG,"%d ",Gf_h[i*dimy+j]);
			fprintf(archB,"%d ",Bf_h[i*dimy+j]);
		}
		fprintf(archTopo,"\n");
		fprintf(archR,"\n");
		fprintf(archG,"\n");
		fprintf(archB,"\n");
	}
	fclose(archTopo);
	fclose(archR);
	fclose(archG);
	fclose(archB);
	*/

	//***************** Archivo de varianza final
	FILE *archVar;
	archVar=fopen("Resultados/VarUltima","w+");
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(archVar,"%f ",var_h[i*dimy+j]);

		}
		fprintf(archVar,"\n");
	}
	fclose(archVar);


	free(var_h);
	free(max_h);
	free(topof_h);
	free(R_h);
	free(G_h);
	free(B_h);
	free(Rf_h);
	free(Gf_h);
	free(Bf_h);

	cudaFree(var_d);
	cudaFree(max_d);
	cudaFree(topof_d);
	cudaFree(R_d);
	cudaFree(G_d);
	cudaFree(B_d);
	cudaFree(Rf_d);
	cudaFree(Gf_d);
	cudaFree(Bf_d);

	
	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	//getchar ();
	return 0;

}//FIN funcion main()
