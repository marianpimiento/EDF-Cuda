/*
The code generates a 3D image of a stack of images.

For each image (matrix) calculate the variance at all points, and then create a topography matrix (relief matrix) with
the position (number in the stack) of the image that had the largest variance in a pixel. The same with the color of the
image (RGB matrices).
*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h>

//************Global variables***************

struct point
{
	int x;
	int y;
};
#define IJ_TO_ID(i,j) (((i)*dimy)+(j))

//************** Kernel CUDA *********************
__global__  void EDF(int *R_d, int *G_d, int *B_d, int *Rf_d, int *Gf_d, int *Bf_d, int *topof_d, long double *max_d, int d) {

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int id = idx + idy*blockDim.x*gridDim.x;

	//int id = idy + idx*blockDim.y*gridDim.y;

	int dimx = 1040, dimy = 1388, tam_imag = 1040*1388, msk = 3, M_d[9], k;
	long double X = 0.f, Xprom = 0.f, Y = 0.f, var = 0.f;
	//Rf_d[id] = id;
	//int img_x = (id) % dimx;
	//int img_y = (id) / dimx;

	int img_x = id / dimy;
	int img_y = id % dimy;

	//int i = 0;
	//unsigned long long int id2;

	M_d[0] = ((img_x < 1 || img_y < 1) ? 0 : G_d[IJ_TO_ID(img_x - 1, img_y - 1)]);
	M_d[1] = ((img_x < 1) ? 0 : G_d[IJ_TO_ID(img_x - 1, img_y)]);
	M_d[2] = ((img_x<1 || img_y>dimy - 2) ? 0 : G_d[IJ_TO_ID(img_x - 1, img_y + 1)]);
	M_d[3] = ((img_y < 1) ? 0 : G_d[IJ_TO_ID(img_x, img_y - 1)]);  //img_x
	M_d[4] = G_d[IJ_TO_ID(img_x, img_y)];
	M_d[5] = ((img_y > dimy - 2) ? 0 : G_d[IJ_TO_ID(img_x, img_y + 1)]);
	M_d[6] = ((img_x > dimx - 2 || img_y < 1) ? 0 : G_d[IJ_TO_ID(img_x + 1, img_y - 1)]);
	M_d[7] = ((img_x > dimx - 2) ? 0 : G_d[IJ_TO_ID(img_x + 1, img_y)]);
	M_d[8] = ((img_x > dimx - 2 || img_y > dimy - 2) ? 0 : G_d[IJ_TO_ID(img_x + 1, img_y + 1)]);


	for (k = 0;k < msk*msk;k++)
		X += M_d[k];

	Xprom = ((long double)X) / (msk*msk);

	for (k = 0;k < msk*msk;k++)
		Y += (Xprom - M_d[k])*(Xprom - M_d[k]);
	
	var = ((long double)Y) / (msk*msk);
	
	//syncthreads();
	__syncthreads();

	if (var > max_d[id]) {
		topof_d[id] = d;
		Rf_d[id] = R_d[id];
		Gf_d[id] = G_d[id];
		Bf_d[id] = B_d[id];
		max_d[id] = var;
	}
}


long msk = 3, dimx = 1040, dimy = 1388, tam_imag = 1040*1388;
//*****************Main function**********************
int main(int argc, char* argv[]) {

	//***************Variables**************
	int i, j, m, cont, tam_B, init, fin;
	cudaError_t cudaStatus;
	FILE *matrizR, *matrizG, *matrizB;
	int d;
	float t;
	clock_t tinicio, t_GPU;
	tinicio = clock();

	int *topof_h, *R_h, *G_h, *B_h, *Rf_h, *Gf_h, *Bf_h;
	long double *max_h;

	int *topof_d, *R_d, *G_d, *B_d, *Rf_d, *Gf_d, *Bf_d;
	long double *max_d;

	//************ Malloc in host and device *************** 
	R_h = (int *)malloc(sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&R_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for R_d!\n");
		exit(0);
	}
	G_h = (int *)malloc(sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&G_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for G_d Line=%d!\n", __LINE__);
		exit(0);
	}
	B_h = (int *)malloc(sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&B_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for R_d!\n");
		exit(0);
	}

	Rf_h = (int *)malloc(sizeof(int)*tam_imag);
	memset((void*)Rf_h, 0, sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&Rf_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for Rf_d!\n");
		exit(0);
	}
	Gf_h = (int *)malloc(sizeof(int)*tam_imag);
	memset((void*)Gf_h, 0, sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&Gf_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed for Gf_d line %d!\n", __LINE__);
		fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
		exit(0);
	}
	Bf_h = (int *)malloc(sizeof(int)*tam_imag);
	memset((void*)Bf_h, 0, sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&Bf_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for Bf_d!\n");
		fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
		exit(0);
	}
	topof_h = (int *)malloc(sizeof(int)*tam_imag);
	memset((void *)topof_h, 0, sizeof(int)*tam_imag);
	cudaStatus = cudaMalloc((void**)&topof_d, tam_imag * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed for topof_d line %d!\n", __LINE__);
		fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
		exit(0);
	}

	//cudaMemset((void *)topof_d, 0, tam_imag * sizeof(int)); //hosam

	max_h = (long double *)malloc(sizeof(long double)*tam_imag);
	memset((void*)max_h, 0, sizeof(long double)*tam_imag);
	cudaStatus = cudaMalloc((void**)&max_d, tam_imag * sizeof(long double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed for max_d line %d!\n", __LINE__);
		fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
		exit(0);
	}
	//cudaMemset(max_d, 0, sizeof(float)*tam_imag);

	//cudaMemset((void *)max_h, 0, sizeof(float)*tam_imag);

	//init=atoi(argv[1]);
	//fin=atoi(argv[2]);
	init = 1;
	fin = 20;

	//*************** Principal FOR ****************
	for (d = init;d <= fin;d++)
	{

		printf("d=%d \n", d);
		//*****************Read RGB files****************
		char rutaR[1024];
		sprintf(rutaR, "%s%d%s", "RGB/", d, "/R");
		matrizR = fopen(rutaR, "r+");
		if (!matrizR)
		{
			printf("Error open file R\n");
			exit(0);
		}

		char rutaG[1024];
		sprintf(rutaG, "%s%d%s", "RGB/", d, "/G");
		matrizG = fopen(rutaG, "r+");
		if (!matrizG)
		{
			printf("Error open file G\n");
			exit(0);
		}

		char rutaB[1024];
		sprintf(rutaB, "%s%d%s", "RGB/", d, "/B");
		matrizB = fopen(rutaB, "r+");
		if (!matrizB)
		{
			printf("Error open file B\n");
			exit(0);
		}

		memset((void*)R_h, 0, sizeof(int)*tam_imag);
		memset((void*)G_h, 0, sizeof(int)*tam_imag);
		memset((void*)B_h, 0, sizeof(int)*tam_imag);

		for(i=0;i<dimx;i++){
			for(j=0;j<dimy;j++){
				fscanf(matrizR, "%d", &R_h[i*dimy + j]);
				fscanf(matrizG, "%d", &G_h[i*dimy + j]);
				fscanf(matrizB, "%d", &B_h[i*dimy + j]);
			}
		}
		fclose(matrizR);
		fclose(matrizG);
		fclose(matrizB);

		//***************** Kernel EDF *******************

		cudaStatus = cudaMemcpy(R_d, R_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for R_d line %d!\n", __LINE__);
			exit(0);
		}

		cudaStatus = cudaMemcpy(G_d, G_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for G_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(B_d, B_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for B_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}

		cudaStatus = cudaMemcpy(Rf_d, Rf_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Rf_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(Gf_d, Gf_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Gf_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(Bf_d, Bf_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Bf_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}

		cudaStatus = cudaMemcpy(topof_d, topof_h, sizeof(int)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for topof_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(max_d, max_h, sizeof(long double)*tam_imag, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for max_d line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}

		dim3 Grid(347, 20);
		dim3 Block(13, 16);

		EDF <<<Grid, Block >>> (R_d, G_d, B_d, Rf_d, Gf_d, Bf_d, topof_d, max_d, d);

		//printf("\n\n FINISH \n\n");
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ The code stops here
		cudaStatus = cudaMemcpy(Rf_h, Rf_d, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Rf_h line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(Gf_h, Gf_d, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Gf_h line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(Bf_h, Bf_d, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for Bf_h line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}

		cudaStatus = cudaMemcpy(topof_h, topof_d, sizeof(int)*tam_imag, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for topof_h line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}
		cudaStatus = cudaMemcpy(max_h, max_d, sizeof(long double)*tam_imag, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for max_h line %d!\n", __LINE__);
			fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(cudaStatus));
			exit(0);
		}

	} //End for

	  //****************Save results**************
	FILE *archTopo, *archR, *archG, *archB;
	archTopo = fopen("Resultados/topo-f3.txt", "w+");
	archR = fopen("Resultados/R-f3.txt", "w+");
	archG = fopen("Resultados/G-f3.txt", "w+");
	archB = fopen("Resultados/B-f3.txt", "w+");
	for(i=0;i<dimx;i++) {
		for(j=0;j<dimy;j++) {
			fprintf(archTopo, "%d ", topof_h[i*dimy + j]);
			fprintf(archR, "%d ", Rf_h[i*dimy + j]);
			fprintf(archG, "%d ", Gf_h[i*dimy + j]);
			fprintf(archB, "%d ", Bf_h[i*dimy + j]);
		}
		fprintf(archTopo, "\n");
		fprintf(archR, "\n");
		fprintf(archG, "\n");
		fprintf(archB, "\n");
	}
	/*for(i=0;i<tam_imag;i++) {

			fprintf(archTopo, "%d ", topof_h[i]);
			fprintf(archR, "%d ", Rf_h[i]);
			fprintf(archG, "%d ", Gf_h[i]);
			fprintf(archB, "%d ", Bf_h[i]);

	}*/

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


	t_GPU = clock();
	t = ((float)t_GPU - (float)tinicio) / CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n", t);

	//getchar ();
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}

	return 0;

}//END Main function
