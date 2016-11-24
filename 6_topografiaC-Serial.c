#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//************variables globales***************
int dimx=1040, dimy=1388, tam_imag=1388*1040, N=20;

//****************declaracion de funciones****************

float** leerMatrizVarianza(int d);

//*****************funcion main**********************

int main(int argc,char* argv[]){

	//***************declaracion de variables**************
	//int N=atoi(argv[2]);
	int i,j,k,temp, **top;
	float **max,**varianza, t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	//************inicializacion de variables***************

	max=(float **)malloc(sizeof(float)*dimx);
	for(i=0;i<dimx;i++)
		max[i]=(float*)malloc(sizeof(float)*dimy);

	top=(int **)malloc(sizeof(int)*dimx);
	for(i=0;i<dimx;i++)
		top[i]=(int*)malloc(sizeof(int)*dimy);

	varianza=(float **)malloc(sizeof(float)*dimx);
	for(i=0;i<dimx;i++)
		varianza[i]=(float*)malloc(sizeof(float)*dimy);

	//***************calculo de la mayor varianza************

	temp=1;
	printf("k=%d\n",temp);
	max=leerMatrizVarianza(temp);
	for(i=0;i<dimx;i++)
		for(j=0;j<dimy;j++)
			top[i][j]=temp;

	for(k=2;k<=N;k++){
		printf("k=%d\n",k);
		varianza=leerMatrizVarianza(k);

		for(i=0;i<dimx;i++){
			for(j=0;j<dimy;j++){
				if(varianza[i][j]>max[i][j]){
					top[i][j]=k;
					max[i][j]=varianza[i][j];
				}
			}
		}
	}

	free(varianza);

	//***************Almacenamiento matriz topografica************

	FILE *topo;
	topo=fopen("Resultados/topo6","w+");
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(topo,"%d ",top[i][j]);
		}
		fprintf(topo,"\n");
	}

	fclose(topo);
	free(max);

	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("tiempo de procesamiento: %6.3f s\n",t);

	getchar ();
	return 0;

}//FIN funcion main()


//******************leerMatrizVarianza****************

float** leerMatrizVarianza(int d){
	int i,j;

	char rutavar[]="VARIANZAS/";
	sprintf(rutavar,"%s%d",rutavar,d);

	FILE* archivo;
	archivo=fopen(rutavar,"r+") ;

	float **var=(float **)malloc(sizeof(float)*dimx);

	for(i=0;i<dimx;i++)
		var[i]=(float*)malloc(sizeof(float)*dimy);

	for(i=0;i<dimx;i++)
		for(j=0;j<dimy;j++){
			fscanf(archivo,"%f",&var[i][j]);
		}

	fclose(archivo);

	return var;
}













