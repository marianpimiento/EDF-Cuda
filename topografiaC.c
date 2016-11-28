//*************inclución de librerias***************

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//************variables globales***************

int N=4, dimx=1388, dimy=1040, tam_imag=1388*1040;

//****************declaración de funciones****************

float** leerMatrizVarianza(int d);

//*****************función main**********************

int main(int argc,char* argv[]){

	//***************declaración de variables**************

	int i,j,k,temp, **top;
	float **max,**varianza, t;
	clock_t tinicio, t_GPU;
	tinicio=clock();

	//************inicialización de variables***************

	max=(float **)malloc(sizeof(float)*dimx);
	for(i=0;i<dimx;i++)
		max[i]=(float*)malloc(sizeof(float)*dimy);

	top=(int **)malloc(sizeof(int)*dimx);
	for(i=0;i<dimx;i++)
		top[i]=(int*)malloc(sizeof(int)*dimy);

	varianza=(float **)malloc(sizeof(float)*dimx);
	for(i=0;i<dimx;i++)
		varianza[i]=(float*)malloc(sizeof(float)*dimy);

	//***************cálculo de la mayor varianza************

	temp=1;
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

	//***************Almacenamiento matriz topográfica************

	FILE *topo;
	topo=fopen("Resultados/topos","w+");
	for(i=0;i<dimx;i++){
		for(j=0;j<dimy;j++){
			fprintf(topo,"%d ",top[i][j]);
			//printf("%d ",top[i][j]);
			if(j%2559==0 && j!=0){
				fprintf(topo,"\n");
				//printf("\n");
			}
		}
	}

	fclose(topo);
	free(max);

	t_GPU=clock();
	t = ((float)t_GPU-(float)tinicio)/CLOCKS_PER_SEC;
	printf("tiempo de procesamiento: %6.3f s\n",t);

	return 0;

}//FIN función main()


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













