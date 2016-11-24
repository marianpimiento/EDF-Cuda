#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//****Variables globales****
int dimx=1040, dimy=1388,tam_imag=1388*1040;

int main(int argc,char* argv[]){

	//****declaracion de variables para el host y device****
	//int N=atoi(argv[2]);
	int N=20;
	int i, j, d, cont;
	int *R_h, *G_h, *B_h, *T_h;
	int *R, *G, *B;

	FILE *file, *Red, *Green, *Blue;
	FILE *ArchivoR, *ArchivoG, *ArchivoB;

	//****Leer archivo de la matriz topografica****

	file=fopen("Resultados/topo6","r+");

	R_h=(int *)malloc(sizeof(int)*tam_imag);
	G_h=(int *)malloc(sizeof(int)*tam_imag);
	B_h=(int *)malloc(sizeof(int)*tam_imag);
	T_h=(int *)malloc(sizeof(int)*tam_imag);

	//******matrices Resultados finales******

	R=(int *)malloc(sizeof(int)*tam_imag);
	G=(int *)malloc(sizeof(int)*tam_imag);
	B=(int *)malloc(sizeof(int)*tam_imag);

	//***********calculo del tiempo de procesamiento***********

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
	fclose(file);

	//******operaciones*******

	for(d=1;d<=N;d++){

		//*******matriz R*********
		printf("d=%d\n", d);
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

		//****************** calculo RGB ****************** +++++++++++++++++++++++++++++++++++++++

		for(i=0;i<tam_imag;i++){
			if(T_h[i]==d){
				R[i]=R_h[i];
				G[i]=G_h[i];
				B[i]=B_h[i];
			}
		}

	}//Fin for

	//almacenamiento de las matrices resultantes*******

	ArchivoR=fopen("Resultados/R7","w+");
	ArchivoG=fopen("Resultados/G7","w+");
	ArchivoB=fopen("Resultados/B7","w+");

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
	printf("tiempo de procesamiento de RGB:%6.3f s\n",t);

	getchar ();
	return 0;
}//Fin funcion main()
