#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<cuda.h>

//************variables globales***************
int dimx=1040, dimy=1388, tam_imag=1388*1040;

//**********KERNEL**************
__global__ void kernel (float *max, float *var, int *top, int k){
	int idx=threadIdx.x + blockIdx.x*blockDim.x;
	int tam_imag=1388*1040;

	if(idx<tam_imag){
		if(var[idx]>max[idx]){
			top[idx]=k;
			max[idx]=var[idx];
		}
	}
}

float *leerMatrizVarianza(int d);

//*****************funcion main**********************

int main(int argc,char* argv[]){

	//***************declaracion de variables**************
	int N=atoi(argv[2]);
	int i,k,temp;
	int *top_d; int top_h[dimx*dimy];
	cudaMalloc((void **)&top_d,sizeof(int)*dimx*dimy);

	float *max_d, *var_d;
	float *max_h, *var_h;

	var_h=(float *)malloc(sizeof(float)*dimx*dimy);
	max_h=(float *)malloc(sizeof(float)*dimx*dimy);
	cudaMalloc((void **)&max_d,sizeof(float)*dimx*dimy);
	cudaMalloc((void **)&var_d,sizeof(float)*dimx*dimy);

	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	//***************calculo de la mayor varianza************

	temp=1;
	max_h=leerMatrizVarianza(temp);
	for(i=0;i<dimx*dimy;i++)
		top_h[i]=temp;

	for(k=2;k<=N;k++){
		printf("k=%d\n", k);
		var_h=leerMatrizVarianza(k);
		cudaMemcpy(max_d,max_h,sizeof(float)*dimx*dimy,cudaMemcpyHostToDevice);
		cudaMemcpy(var_d,var_h,sizeof(float)*dimx*dimy,cudaMemcpyHostToDevice);
		cudaMemcpy(top_d,top_h,sizeof(int)*dimx*dimy,cudaMemcpyHostToDevice);

		kernel<<<6940,208>>>(max_d,var_d,top_d,k);

		cudaMemcpy(top_h,top_d,sizeof(int)*dimx*dimy,cudaMemcpyDeviceToHost);
		cudaMemcpy(max_h,max_d,sizeof(float)*dimx*dimy,cudaMemcpyDeviceToHost);
	}
	
	
	cudaFree(max_d);
	cudaFree(var_d);
	cudaFree(top_d);

	FILE *topo;
	topo=fopen("Resultados/topo","w+");
	for(i=0;i<dimx*dimy;i++){
		if(i%dimy==0 && i!=0)
			fprintf(topo,"\n");
			fprintf(topo,"%d ",top_h[i]);
	}

	fclose(topo);
	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("tiempo de procesamiento de topografia: %6.3f s\n",t);

}//FIN funcion main()


//******************leerMatrizVarianza****************

float* leerMatrizVarianza(int d){
	int i;
	char rutavar[]="VARIANZAS/";
	sprintf(rutavar,"%s%d",rutavar,d);

	FILE* archivo;
	archivo=fopen(rutavar,"r") ;

	float *var;
	var=(float *)malloc(sizeof(float)*dimx*dimy);

	for(i=0;i<dimx*dimy;i++)
		fscanf(archivo,"%f",&var[i]);
	fclose(archivo);

	return var;

}