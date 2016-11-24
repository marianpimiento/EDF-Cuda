#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<time.h>

//************variables globales***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;

//*******************kernel********************

__global__ void kernel (int *Gext_d,float *var_d){

	int i, dimy_ext, id_p, M_d[9], dimy=1388,tam_imag=1388*1040,msk=3;
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	int idy = threadIdx.y + blockIdx.y*blockDim.y;
	int offset=idx + idy*blockDim.x*gridDim.x;
	int id=offset;

	float X=0.f,Xprom=0.f,Y=0.f;
		
	var_d[id]=0;
	dimy_ext=dimy+2;

	if(offset<tam_imag){
		
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
		var_d[id]=Y/(msk*msk);

	}
}


//*****************Funcion Main**********************

int main(int argc,char* argv[]){

	//***************Declaracion de variables**************

	int i,j,d,m,cont,tam_ext,init,fin;
	int *Gext_d, *Gext_h;
	//, **G;

	float t, *var_d,*var_h;

	FILE *archV, *matrizG;

	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	clock_t tinicio, t_GPU;
	tinicio=clock();

	tam_ext=(dimx+2)*(dimy+2);

	//*************For que recorre todas las imagenes ************
	for(d=init;d<=fin;d++){

		printf("d=%d \n", d);
		
		Gext_h=(int *)malloc(sizeof(int)*tam_ext);
		cudaMalloc((void**)&Gext_d, tam_ext*sizeof(int));
		var_h=(float *)malloc(sizeof(float)*tam_imag);
		cudaMalloc((void**)&var_d,tam_imag*sizeof(float));
		

		//*****************Lectura Matriz G****************

		char rutaG[]="MiTesis/";
		sprintf(rutaG, "%s%d%s","RGB/",d,"/G"); 
		matrizG=fopen(rutaG,"r+");

		cont=0;
		for(i=0;i<dimx+2;i++){
			for(j=0;j<dimy+2;j++){
				if (i==0 || j==0 || i==dimx+1 || j==dimy+1){
					Gext_h[cont]=0;
				} else{
					fscanf(matrizG, "%d", &Gext_h[cont]); 
				}
				cont++;
			}
		}
		fclose(matrizG);


		//******************Llamado de kernel*******************

		dim3 Grid(347,20);
		dim3 Block(13,16);

		cudaMemcpy(Gext_d,Gext_h,sizeof(int)*tam_ext,cudaMemcpyHostToDevice);

		kernel<<<Grid,Block>>>(Gext_d,var_d);

		cudaMemcpy(var_h,var_d,sizeof(float)*tam_imag,cudaMemcpyDeviceToHost);


		//****************Almacenamiento matriz de Varianza**************

		char rutaV[]="VARIANZAS/";
		sprintf(rutaV, "%s%d", rutaV,d);
		archV=fopen(rutaV,"w+");

		for(m=0;m<tam_imag;m++){
			if(m%dimy==0 && m!=0){
				fprintf(archV,"\n");
			}
			fprintf(archV,"%f ",var_h[m]);
		}
		fclose(archV);

		free(Gext_h);
		free(var_h);
		//free(G);
		cudaFree(var_d);
		cudaFree(Gext_d);
	}


	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	return 0;

}//FIN funcion main()
