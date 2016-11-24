#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<time.h>

//************variables globales***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;

//*******************kernel********************

__global__ void varianza (int *Gext_d,float *var_d){

	int i, dimy_ext, id_p, M_d[9], dimy=1388,tam_imag=1388*1040,msk=3;
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int offset=idx + idy*blockDim.x*gridDim.x;
	int id=offset;

	float X=0.f,Xprom=0.f,Y=0.f;
		
	//float var=0;
	//var_d[id]=0;
	
	if(offset<tam_imag){
		
		dimy_ext=dimy+2;
		offset+=2*idy;
		id_p=offset+(dimy+msk);

		M_d[0]=Gext_d[offset];
		M_d[1]=Gext_d[offset+1];
		M_d[2]=Gext_d[offset+2];
		M_d[3]=Gext_d[id_p-1];
		M_d[4]=Gext_d[id_p];
		M_d[5]=Gext_d[id_p+1];
		M_d[6]=Gext_d[(id_p-1)+dimy_ext];
		M_d[7]=Gext_d[id_p+dimy_ext];
		M_d[8]=Gext_d[(id_p+1)+dimy_ext];

		for(i=0;i<msk*msk;i++)
			X+=M_d[i];
		Xprom=((float)X)/(msk*msk);

		for(i=0;i<msk*msk;i++)
			Y+=(Xprom-M_d[i])*(Xprom-M_d[i]);
		
		//var=Y/(msk*msk);
		var_d[id]=Y/(msk*msk);

	}
}


__global__ void topografia (float *var_d,int *topof_d,float *max_d, int d){

	int idx=threadIdx.x + blockIdx.x*blockDim.x;
	int tam_imag=1388*1040;

	if(idx<tam_imag){
		if(var_d[idx]>max_d[idx]){
			topof_d[idx]=d;
			max_d[idx]=var_d[idx];
			/*	Rf_d[id]=R_d[id];
			Gf_d[id]=G_d[id];
			Bf_d[id]=B_d[id];*/
		}
	}
}

//*****************Funcion Main**********************

int main(int argc,char* argv[]){

	//***************Declaracion de variables**************

	int i,j,d,m,cont,tam_ext,init,fin;

	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	FILE *matrizR, *matrizG, *matrizB, *matrizGext;
	
	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	tam_ext=(dimx+2)*(dimy+2);

	int *topof_h, *R_h, *G_h, *B_h, *Rf_h, *Gf_h, *Bf_h, *Gext_h;
	float *max_h, *var_h;

	int *topof_d, *R_d, *G_d, *B_d, *Rf_d, *Gf_d, *Bf_d, *Gext_d;
	float *max_d, *var_d;


	//************Inicializacion de variables en el host y en el device *************** 
	
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

	Gext_h=(int *)malloc(sizeof(int)*tam_ext);
	cudaMalloc((void**)&Gext_d, tam_ext*sizeof(int));

	max_h=(float *)malloc(sizeof(float)*tam_imag);
	cudaMalloc((void**)&max_d, tam_imag*sizeof(float));
	
	//cudaMemset((void *) max_d, 0, sizeof(float)*tam_imag);
	//void *memset(void *str, int c, size_t n)
	//memset((void *) max_h, 0, sizeof(float)*tam_imag);

	for(i=0;i<tam_imag;i++){
		max_h[i]=0.0;
		topof_h[i]=0;
	}

	printf("Antes for principal\n");

	//*************For que recorre todas las imagenes ************
	for(d=init;d<=fin;d++){

		printf("d=%d \n", d);

		var_h=(float *)malloc(sizeof(float)*tam_imag);
		cudaMalloc((void**)&var_d,tam_imag*sizeof(float));

		for(i=0;i<tam_imag;i++){
			var_h[i]=0;
		}
		
		//*****************Lecura de matrices RGB en el host****************
/*
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
*/
		//G extendido

		char rutaGext[]="";
		sprintf(rutaGext, "%s%d%s","RGB/",d,"/G"); 
		matrizGext=fopen(rutaGext,"r+");

		cont=0;
		for(i=0;i<dimx+2;i++){
			for(j=0;j<dimy+2;j++){
				if (i==0 || j==0 || i==dimx+1 || j==dimy+1){
					Gext_h[cont]=0;
				} else{
					fscanf(matrizGext, "%d", &Gext_h[cont]); 
				}
				cont++;
			}
		}
		fclose(matrizGext);

		printf("Despues lectura matrices \n");


		//******************Llamado kernel varianza*******************  ++++++++++++++++++++++++++++++++++++
		printf("*Kenel varianza \n");
		cudaMemcpy(Gext_d,Gext_h,sizeof(int)*tam_ext,cudaMemcpyHostToDevice);

		printf("Despues copia a device\n");

		dim3 Grid(347,20);
		dim3 Block(13,16);

		varianza<<<Grid,Block>>>(Gext_d,var_d);
		printf("Despues kernel \n");

		cudaMemcpy(var_h,var_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);
		printf("Despues copia a host\n");
		printf("var_h[0]= %f\n", var_h[0]);


		//******************Llamado kernel topografia******************* ++++++++++++++++++++++++++++++++
		printf("*Kenel topografia \n");
/*
		cudaMemcpy(R_d,R_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(G_d,G_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(B_d,B_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);

		cudaMemcpy(Rf_d,Rf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Gf_d,Gf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(Bf_d,Bf_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
*/
		cudaMemcpy(var_d,var_h,sizeof(float)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(topof_d,topof_h,sizeof(int)*tam_imag,cudaMemcpyHostToDevice);
		cudaMemcpy(max_d,max_h,sizeof(float)*tam_imag,cudaMemcpyHostToDevice);

		printf("Despues copia a device\n");

		//dim3 Grid(347,20);
		//dim3 Block(13,16);

		//topografia<<<6940,208>>>(R_d,G_d,B_d,Rf_d,Gf_d,Bf_d,topof_d,max_d,var_d,d);
		topografia<<<6940,208>>>(var_d,topof_d,max_d,d);
		printf("Despues kernel \n");

/*		cudaMemcpy(Rf_h,Rf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Gf_h,Gf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(Bf_h,Bf_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
*/
		cudaMemcpy(topof_h,topof_d,sizeof(int)*tam_imag,cudaMemcpyDeviceToHost);
		cudaMemcpy(max_h,max_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);
		printf("Despues copia a host\n");
		printf("topof_h[0]= %d\n", topof_h[0]);

	}//Finaliza For principal

	//****************Almacenamiento matrices**************

	FILE *archTopo, *archR, *archG, *archB, *archV;
	archTopo=fopen("Resultados/topos12","w+");
	/*archR=fopen("Resultados/R12","w+");
	archG=fopen("Resultados/G12","w+");
	archB=fopen("Resultados/B12","w+");
	archV=fopen("Resultados/VarUltima","w+");*/
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(archTopo,"%d ",topof_h[i*dimy+j]);
			/*fprintf(archR,"%d ",Rf_h[i*dimy+j]);
			fprintf(archG,"%d ",Gf_h[i*dimy+j]);
			fprintf(archB,"%d ",Bf_h[i*dimy+j]);
			fprintf(archV,"%f ",var_h[i*dimy+j]);*/
		}
		fprintf(archTopo,"\n");
		/*fprintf(archR,"\n");
		fprintf(archG,"\n");
		fprintf(archB,"\n");
		fprintf(archV,"\n");*/
	}
	fclose(archTopo);
/*	fclose(archR);
	fclose(archG);
	fclose(archB);
	fclose(archV);*/


	//****************Libera memoria**************
	free(R_h);
	cudaFree(R_d);
	free(G_h);
	cudaFree(G_d);
	free(B_h);
	cudaFree(B_d);
	free(Rf_h);
	cudaFree(Rf_d);
	free(Gf_h);
	cudaFree(Gf_d);
	free(Bf_h);
	cudaFree(Bf_d);
	free(Gext_h);
	cudaFree(Gext_d);
	free(topof_h);
	cudaFree(topof_d);
	free(max_h);
	cudaFree(max_d);

	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	return 0;

}//FIN funcion main()
