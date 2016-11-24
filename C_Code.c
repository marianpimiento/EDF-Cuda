/*
The code generates a 3D image of a stack of images.

For each image (matrix) calculate the variance at all points, and then create a topography matrix (relief matrix) with 
the position (number in the stack) of the image that had the largest variance in a pixel. The same with the color of the 
image (RGB matrices).
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//************Global variables***************
int msk=3, dimx=1040, dimy=1388, tam_imag=1388*1040;

//*****************Main function*********************
int main(int argc,char* argv[]){

	//***************Viariables**************
	int i,j,m,cont,tam_B, init,fin;
	init=atoi(argv[1]);
	fin=atoi(argv[2]);

	FILE *arch, *matrizR, *matrizG, *matrizB, *archM;
	int d;
	float t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	int **topof, **R, **G, **B, **Rf, **Gf, **Bf;
	float **max, *var_h;
	
	
	max=(float **)malloc(sizeof(float)*dimx);
	topof=(int **)malloc(sizeof(int)*dimx);
	R=(int **)malloc(sizeof(int)*dimx);
	G=(int **)malloc(sizeof(int)*dimx);
	B=(int **)malloc(sizeof(int)*dimx);
	Rf=(int **)malloc(sizeof(int)*dimx);
	Gf=(int **)malloc(sizeof(int)*dimx);
	Bf=(int **)malloc(sizeof(int)*dimx);
	for(i=0;i<dimx;i++){
		max[i]=(float*)malloc(sizeof(float)*dimy);
		topof[i]=(int*)malloc(sizeof(int)*dimy);
		R[i]=(int*)malloc(sizeof(int)*dimy);
		G[i]=(int*)malloc(sizeof(int)*dimy);
		B[i]=(int*)malloc(sizeof(int)*dimy);
		Rf[i]=(int*)malloc(sizeof(int)*dimy);
		Gf[i]=(int*)malloc(sizeof(int)*dimy);
		Bf[i]=(int*)malloc(sizeof(int)*dimy);
	}
	
	var_h=(float *)malloc(sizeof(float)*tam_imag);
	
	//***************Principal FOR****************
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
				fscanf(matrizR, "%d", &R[i][j]);
				fscanf(matrizG, "%d", &G[i][j]); 
				fscanf(matrizB, "%d", &B[i][j]); 
			}
		}
		fclose(matrizR);
		fclose(matrizG);
		fclose(matrizB);

		//********************* Variance ********************
		
		int M_d[9], k;
		float X,Xprom,Y;
		
		for(i=0;i<dimx;i++){
			for(j=0; j<dimy; j++){
				X=0.f;
				Xprom=0.f;
				Y=0.f;
				
				var_h[i*dimy+j]=0;
				M_d[0]=((i<1 || j<1) ? 0:G[i-1][j-1]);
				M_d[1]=((i<1) ? 0:G[i-1][j]);
				M_d[2]=((i<1 || j>dimy-2) ? 0:G[i-1][j+1]);
				M_d[3]=((j<1) ? 0:G[i][j-1]);
				M_d[4]=G[i][j];
				M_d[5]=((j>dimy-2) ? 0:G[i][j+1]);
				M_d[6]=((i>dimx-2 || j<1) ? 0:G[i+1][j-1]);
				M_d[7]=((i>dimx-2) ? 0:G[i+1][j]);
				M_d[8]=((i>dimx-2 || j>dimy-2) ? 0:G[i+1][j+1]);

				
				for(k=0;k<msk*msk;k++){
					X+=M_d[k];
				}
				
				Xprom=((float)X)/(msk*msk);

				for(k=0;k<msk*msk;k++)
					Y+=(Xprom-M_d[k])*(Xprom-M_d[k]);

				var_h[i*dimy+j]=Y/(msk*msk);

			}
		}


		//**********Topography********************* 

		if(d==init){
			for(i=0;i<dimx;i++){
				for(j=0;j<dimy;j++){
					max[i][j]=var_h[i*dimy+j];
					topof[i][j]=init;
					Rf[i][j]=R[i][j];
					Gf[i][j]=G[i][j];
					Bf[i][j]=B[i][j];
				}
			}
		}
		
		if(d>init){
			for(i=0;i<dimx;i++){
				for(j=0;j<dimy;j++){
					if(var_h[i*dimy+j]>max[i][j]){
						topof[i][j]=d;
						Rf[i][j]=R[i][j];
						Gf[i][j]=G[i][j];
						Bf[i][j]=B[i][j];
						max[i][j]=var_h[i*dimy+j];
					}
				}
			}
		}

	}

	FILE *archTopo, *archR, *archG, *archB;
	archTopo=fopen("Resultados/topos9","w+");
	archR=fopen("Resultados/R9","w+");
	archG=fopen("Resultados/G9","w+");
	archB=fopen("Resultados/B9","w+");
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(archTopo,"%d ",topof[i][j]);
			fprintf(archR,"%d ",Rf[i][j]);
			fprintf(archG,"%d ",Gf[i][j]);
			fprintf(archB,"%d ",Bf[i][j]);
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

	free(var_h);
	free(max);
	free(topof);
	free(R);
	free(G);
	free(B);
	free(Rf);
	free(Gf);
	free(Bf);

	
	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("\ntiempo de procesamiento de varianzas: %6.3fs\n",t);

	//getchar ();
	return 0;

}//FIN funcion main()
