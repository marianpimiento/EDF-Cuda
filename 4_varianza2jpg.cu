#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<time.h>

//************variables globales***************
int msk=3, dimx=1388, dimy=1040, tam_imag=1388*1040;

//*******************kernel********************

__global__ void kernel (int *B_d,float *var_d){

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int offset=idx + idy*blockDim.x*gridDim.x;

	int id=offset;
	int i;
	float X=0.f,Xprom=0.f,Y=0.f;
	int dimy=1040,tam_imag=1388*1040,msk=3;
	var_d[id]=0;

	if(offset<tam_imag){
		int dimy_B=dimy+2;

		offset+=2*idy;
		int id_p=offset+(dimy+msk);

		int M_d[9];

		M_d[0]=B_d[offset];
		M_d[1]=B_d[offset+1];
		M_d[2]=B_d[offset+2];
		M_d[3]=B_d[id_p-1];
		M_d[4]=B_d[id_p];
		M_d[5]=B_d[id_p+1];
		M_d[6]=B_d[(id_p-1)+dimy_B];
		M_d[7]=B_d[id_p+dimy_B];
		M_d[8]=B_d[(id_p+1)+dimy_B];

		for(i=0;i<msk*msk;i++)
			X+=M_d[i];
		Xprom=((float)X)/(msk*msk);

		for(i=0;i<msk*msk;i++)
			Y+=(Xprom-M_d[i])*(Xprom-M_d[i]);
		var_d[id]=Y/(msk*msk);

	}
}


//*****************funcion main**********************

int main(int argc,char* argv[]){

	//***************declaracion de variables**************

	printf("\nInicia main");
	int i,j,m,cont,tam_B, init,fin;
	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	tam_B=(dimx+2)*(dimy+2);

	printf("\nCrea FILE\n");
	FILE *arch, *matrizG, *archM; //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	int **A;
	int B[dimx+2][dimy+2];

	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	int *B_d, *B_h;
	float *var_d,*var_h;

	int d;
	for(d=init;d<=fin;d++){

		//*******************declaracion de variables***************
		printf("d=%d \n", d);
		
		B_h=(int *)malloc(sizeof(int)*tam_B);
		cudaMalloc((void**)&B_d, tam_B*sizeof(int));
		var_h=(float *)malloc(sizeof(float)*tam_imag);
		cudaMalloc((void**)&var_d,tam_imag*sizeof(float));
		
		A=(int **)malloc(sizeof(int)*dimx);
		for(i=0;i<dimx;i++)
			A[i]=(int*)malloc(sizeof(int)*dimy);

		//*****************calculo matriz B****************

		char ruta1[]="MiTesis/";
		sprintf(ruta1, "%s%d%s","RGB/",d,"/G2.jpg");  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		matrizG=fopen(ruta1,"r+");

		for(i=0;i<dimx;i++){
			for(j=0;j<dimy;j++){
				fscanf(matrizG, "%d", &A[i][j]); 
			}
		}
		fclose(matrizG);

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		char rutaM[]="matrizP/";
		sprintf(rutaM, "%s%d", rutaM,d);
		archM=fopen(rutaM,"w+");

		for(i=0;i<dimx;i++){
			fprintf(archM,"\n");
			for(j=0;j<dimy;j++){
				fprintf(archM,"%d ",A[i][j]);
			}
		}
		fclose(archM);

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		cont=0;
		for(i=0;i<dimx+2;i++){
			//printf("\n");
			for(j=0;j<dimy+2;j++){
				B[i][j]=((i==0 || j==0 || i==dimx+1 || j==dimy+1) ? 0:A[i-1][j-1]);
				B_h[cont]=B[i][j];
				cont++;
			}
		}

		//******************llamado de kernel*******************

		dim3 Grid(347,20);
		dim3 Block(13,16);

		cudaMemcpy(B_d,B_h,sizeof(int)*tam_B,cudaMemcpyHostToDevice);

		kernel<<<Grid,Block>>>(B_d,var_d);

		cudaMemcpy(var_h,var_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);


		//****************almacenamiento matriz de varianza**************

		char rutaV[]="VARIANZAS/";
		sprintf(rutaV, "%s%dj", rutaV,d); //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		arch=fopen(rutaV,"w+");

		for(m=0;m<tam_imag;m++){
			if(m%dimy==0 && m!=0){
				fprintf(arch,"\n");
			}
			fprintf(arch,"%f ",var_h[m]); // "%.2f " - Imprimiria 2 decimales
		}
		fclose(arch);

		free(B_h);
		free(var_h);
		free(A);
		cudaFree(var_d);
		cudaFree(B_d);
	}


	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	return 0;

}//FIN funcion main()
