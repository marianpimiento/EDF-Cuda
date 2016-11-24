/*
The code generates a 3D image of a stack of images.

For each image (matrix) calculate the variance at all points, and then create a topography matrix (relief matrix) with 
the position (number in the stack) of the image that had the largest variance in a pixel. The same with the color of the 
image (RGB matrices).
*/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h>

//************Global variables***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;


//************** Kernel CUDA *********************
__global__ void EDF(int *R_d,int *G_d,int *B_d,int *Rf_d,int *Gf_d,int *Bf_d,int *topof_d,float *max_d, int d){

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int id = idx + idy*blockDim.x*gridDim.x;

	int dimx=1040, dimy=1388, tam_imag=1388*1040, msk=3, M_d[9], k;
	float X=0.f,Xprom=0.f,Y=0.f, var=0.f;

	if(id<tam_imag){

		M_d[0]=((idx==0 || idy==0) ? 0:G_d[id-1-dimy]);
		M_d[1]=((idx==0) ? 0:G_d[id-dimy]);
		M_d[2]=((idx==0 || idy==dimy) ? 0:G_d[id+1-dimy]);
		M_d[3]=((idy==0) ? 0:G_d[id-1]);
		M_d[4]=G_d[id];
		M_d[5]=((idy==dimy) ? 0:G_d[id+1]);
		M_d[6]=((idx==dimx || idy==0) ? 0:G_d[id-1+dimy]);
		M_d[7]=((idx==dimx) ? 0:G_d[id+dimy]);
		M_d[8]=((idx==dimx || idy==dimy) ? 0:G_d[id+1+dimy]);

		for(k=0;k<msk*msk;k++)
			X+=M_d[k];
				
		Xprom=((float)X)/(msk*msk);

		for(k=0;k<msk*msk;k++)
			Y+=(Xprom-M_d[k])*(Xprom-M_d[k]);

		var=Y/(msk*msk);

		syncthreads();

		if(var>max_d[id]){
			topof_d[id]=d;
			Rf_d[id]=R_d[id];
			Gf_d[id]=G_d[id];
			Bf_d[id]=B_d[id];
			max_d[id]=var;
		}
	}
}



//*****************Main function**********************
int main(int argc,char* argv[]){

	//***************Variables**************
	int i,j,m,cont,tam_B, init,fin;
	
	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	FILE *matrizR, *matrizG, *matrizB;
	int d;
	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	int *topof_h, *R_h, *G_h, *B_h, *Rf_h, *Gf_h, *Bf_h;
	float *max_h;

	int *topof_d, *R_d, *G_d, *B_d, *Rf_d, *Gf_d, *Bf_d;
	float *max_d;
	
	//************ Malloc in host and device *************** 
	
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

	cudaMemset((void *) max_h, 0, sizeof(float)*tam_imag);
	
	//*************** Principal FOR ****************
	for(d=init;d<=fin;d++){

		printf("d=%d \n", d);
		
		//*****************Read RGB files****************

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

		//***************** Kernel EDF *******************

		cudaMemcpy(R_d,R_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(G_d,G_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(B_d,B_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		cudaMemcpy(Rf_d,Rf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Gf_d,Gf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Bf_d,Bf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		cudaMemcpy(topof_d,topof_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(max_d,max_h,sizeof(float)*tam_imag,cudaMemcpyHostToDevice);

		dim3 Grid(347,20);
		dim3 Block(13,16);

		EDF<<<Grid,Block>>>(R_d,G_d,B_d,Rf_d,Gf_d,Bf_d,topof_d,max_d,d);

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ The code stops here
		cudaMemcpy(Rf_h,Rf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Gf_h,Gf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Bf_h,Bf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);

		cudaMemcpy(topof_h,topof_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(max_h,max_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);

	} //End for

	//****************Save results**************

	FILE *archTopo, *archR, *archG, *archB;
	archTopo=fopen("Resultados/topo","w+");
	archR=fopen("Resultados/R","w+");
	archG=fopen("Resultados/G","w+");
	archB=fopen("Resultados/B","w+");
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

	free(max_h);
	free(topof_h);
	free(R_h);
	free(G_h);
	free(B_h);
	free(Rf_h);
	free(Gf_h);
	free(Bf_h);

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

}//END Main function
