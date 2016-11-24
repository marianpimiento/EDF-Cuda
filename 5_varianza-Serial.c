
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//************variables globales***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;


//*****************funcion main**********************

int main(int argc,char* argv[]){

	//***************declaracion de variables**************

	int i,j,m,cont,tam_B, init,fin;
	//init=atoi(argv[1]);
	//fin=atoi(argv[2]);

	init=1; //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	fin=328;

	FILE *arch, *matrizG, *archM;

	int **A;

	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	float *var_h;

	int d;
	for(d=init;d<=fin;d++){

		//*******************declaracion de variables***************
		printf("d=%d \n", d);
		
		var_h=(float *)malloc(sizeof(float)*tam_imag);
		
		A=(int **)malloc(sizeof(int)*dimx);
		for(i=0;i<dimx;i++)
			A[i]=(int*)malloc(sizeof(int)*dimy);

		//*****************calculo matriz B****************

		char ruta1[]="MiTesis/";
		sprintf(ruta1, "%s%d%s","RGB/",d,"/G"); 
		matrizG=fopen(ruta1,"r+");

		for(i=0;i<dimx;i++){
			for(j=0;j<dimy;j++){
				fscanf(matrizG, "%d", &A[i][j]); 
			}
		}
		fclose(matrizG);

		//*********************Calculo de la Varianza ******************** ++++++++++++++++++++++++++++++
		
		int M_d[9], k;
		float X,Xprom,Y;

		for(i=0;i<dimx;i++){
			for(j=0; j<dimy; j++){
				X=0.f;
				Xprom=0.f;
				Y=0.f;
				
				var_h[i*dimy+j]=0;
				M_d[0]=((i<1 || j<1) ? 0:A[i-1][j-1]);
				M_d[1]=((i<1) ? 0:A[i-1][j]);
				M_d[2]=((i<1 || j>dimy-2) ? 0:A[i-1][j+1]);
				M_d[3]=((j<1) ? 0:A[i][j-1]);
				M_d[4]=A[i][j];
				M_d[5]=((j>dimy-2) ? 0:A[i][j+1]);
				M_d[6]=((i>dimx-2 || j<1) ? 0:A[i+1][j-1]);
				M_d[7]=((i>dimx-2) ? 0:A[i+1][j]);
				M_d[8]=((i>dimx-2 || j>dimy-2) ? 0:A[i+1][j+1]);

				
				for(k=0;k<msk*msk;k++){
					X+=M_d[k];
				}
				

				Xprom=((float)X)/(msk*msk);

				for(k=0;k<msk*msk;k++)
					Y+=(Xprom-M_d[k])*(Xprom-M_d[k]);

				var_h[i*dimy+j]=Y/(msk*msk);

			}
		}

		//****************almacenamiento matriz de varianza**************

		char rutaV[]="VARIANZAS/";
		sprintf(rutaV, "%s%d", rutaV,d);
		arch=fopen(rutaV,"w+");

		for(m=0;m<tam_imag;m++){
			if(m%dimy==0 && m!=0){
				fprintf(arch,"\n");
			}
			fprintf(arch,"%f ",var_h[m]);
		}
		fclose(arch);

		free(var_h);
		free(A);
	}


	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	getchar ();
	return 0;

}//FIN funcion main()
