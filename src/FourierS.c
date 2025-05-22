#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
/*
    IMPORTANTE
Las muestras que se les da a las funciones son los valores x[n] de la DFT y x(t) de la CFT, si no por ejemplo
en la CFT necesitariamos definir una funcion f que tome valores en el tiempo para poder hacer el CFT, para más
simplicidad, he supuesto que estos valores son los que tienen mis muestras, ademas de acotarlas entre [-1,1].

*/

#define PI 3.14159265358979323846
//Las muestras están entre -1 y 1
#define T_MAX 1
#define T_MIN -1



//Fourier discreto 😇😇
void DFT(double complex *Fourier,const double *muestras, const int n){
    for (int i=0;i<n;i++){
        Fourier[i] = 0.0 + 0.0*I;
        for (int j=0;j<n;j++){
            double angle = -2.0*PI*i*j/n;
            Fourier[i] += muestras[j]*cexp(I*angle);
        }
    }    
}

//Fourier continuo 🤯🤯 Suma de rectangulos
void CFT(double complex *Fourier, const double *muestras, const int n, const double paso_temporal){
    for (int i=0;i<n;i++){
        Fourier[i] = 0.0 + 0.0*I;
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega
        //Aqui ya entra en juego tanto el intervalo del tiempo como los valores que dan la funcion, es decir el tiempo esta entre -1 y 1
        //y los valores de la funcion están en mis muestras
        //Vamos desde el minimo hasta el maximo pero con nuestro paso temporal para tomar fourier lo más preciso posible
        for (double j = T_MIN;j<T_MAX; j= j+paso_temporal){
            //printf("Estoy en el segundo FOR");
            //Dudas a como acceder con el paso temporal, preguntar a Pablo
            int indice = (int)((j - T_MIN)/paso_temporal);
            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*j) * paso_temporal;
            }else if (indice >= n-1){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice]*cexp(-I*omega*j) * paso_temporal;
            }
            

        }
    }    
}

//Integral con simpson (2)

void CFT_Simpson(double complex *Fourier, const double *muestras, const int n, const double paso_temporal) {
    for (int i = 0; i < n; i++) {
        Fourier[i] = 0.0 + 0.0 * I;
        double omega = 2.0 * PI * i / (T_MAX - T_MIN);

        for (double j = T_MIN; j < T_MAX - paso_temporal; j += paso_temporal) {
            int indice1 = (int)((j - T_MIN) / paso_temporal);
            int indice2 = indice1 + 1;

            if (indice2 >= n) indice2 = n - 1;

            double x_medio = j + paso_temporal / 2.0;
            int indice_medio = (int)((x_medio - T_MIN) / paso_temporal);
            if (indice_medio >= n) indice_medio = n - 1;

            double simpson = (muestras[indice1] + 4.0 * muestras[indice_medio] + muestras[indice2]) / 6.0;

            Fourier[i] += simpson * cexp(-I * omega * j) * paso_temporal;
        }
    }
}


//Integral con trapecio (3)
void CFT_Trapecio(double complex *Fourier, const double *muestras, const int n, const double paso_temporal) {
    for (int i = 0; i < n; i++) {
        Fourier[i] = 0.0 + 0.0 * I;
        double omega = 2.0 * PI * i / (T_MAX - T_MIN);

        for (double j = T_MIN; j < T_MAX - paso_temporal; j += paso_temporal) {
            int indice1 = (int)((j - T_MIN) / paso_temporal);
            int indice2 = indice1 + 1;

            if (indice2 >= n) indice2 = n - 1; 

            double promedio = (muestras[indice1] + muestras[indice2]) / 2.0;
            Fourier[i] += promedio * cexp(-I * omega * j) * paso_temporal;
        }
    }
}


//Apartador de archivos de entrada y saluda
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char * MILLON = "txt/MuestraGenerada.txt"; //una muestra de un millo de elementos 🫨🫨
const char * TREINTAMIL = "txt/MuestraGenerada30000.txt"; //una muestra de un 30.000 elementos 🫨🫨
const char * CINCUENTAMIL = "txt/MuestraGenerada50000.txt";
const char * CIENMIL = "txt/MuestraGenerada100000.txt";
const char * CIENTOCINCUENTAMIL = "txt/MuestraGenerada150000.txt";
const char * DOSCIENTOSCINCUENTAMIL = "txt/MuestraGenerada250000.txt";
const char * QUINIENTOSMIL = "txt/MuestraGenerada500000.txt";
const char * DISTINTAS20 = "txt/muestras.txt"; //20.000 muestras de elementos aleatorios 😎😎
const char * FUNCIONA = "txt/funciona.txt"; //muestras de funcionamiento 😎😎
const char * FUNCIONA2 = "txt/funciona2.txt"; //muestras de funcionamiento 😎😎
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char * SALIDA = "Resultados/SecuencialDFT.txt"; //salida secuencial (cualquier caso) 😇😇
const char * SALIDA2 = "txt/SecuencialDFT2.txt"; //salida secuencial para el caso de 20.000 muestras 😇😇
const char * SALIDA_CONTINUO = "Resultados/SecuencialCFTNormal.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "Resultados/SecuencialCFTSimpson.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO3 = "Resultados/SecuencialCFTTrapecio.txt"; //salida continuo 🤯🤯
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (){
    FILE * entrada = fopen (QUINIENTOSMIL,"r");
    FILE * salida = fopen(SALIDA,"w");
    FILE * salida_cont = fopen(SALIDA_CONTINUO,"w");
    FILE * salida_cont2 = fopen(SALIDA_CONTINUO2,"w");
    FILE * salida_cont3 = fopen(SALIDA_CONTINUO3,"w");

    if (!entrada){
        printf("Error: No se pudo abrir el archivo de entrada\n");
        exit(1);
    }
    if (!salida){
        printf("Error: No se pudo abrir el archivo de salida\n");
        exit(1);
    }

    if (!salida_cont){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    if (!salida_cont2){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    if (!salida_cont3){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    
    int Tam_Vector_muestras;
    while (fscanf(entrada,"%d",&Tam_Vector_muestras) == 1){
        
        double *muestras = malloc(Tam_Vector_muestras*sizeof(double));
        double complex *Fourier = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont2 = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont3 = malloc(Tam_Vector_muestras*sizeof(double complex));
        
        if (!muestras){
            printf("Error: No se pudo asignar memoria muestras\n");
            exit(1);
        }
        if (!Fourier){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont2){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont3){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }

        for (int i=0;i<Tam_Vector_muestras;i++){
            if (fscanf(entrada,"%lf",&muestras[i]) != 1){
                printf("Error: archivo de entrada de la muestras %d\n",i);
                exit(1);
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de DFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        clock_t inicio = clock();
        DFT(Fourier,muestras,Tam_Vector_muestras);
        clock_t fin = clock();
        double tiempo = (double)(fin-inicio)/CLOCKS_PER_SEC;
        //El tiempo se pasa en ms
        fprintf(salida,"%d %lf\n",Tam_Vector_muestras,tiempo);
        //Para mostrar el vector, esta dentro del archivo 😴😴
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida,"%lf %lf\n",creal(Fourier[i]),cimag(Fourier[i]));
        }*/

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de CFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //printf("valor tammuestras %d\n",Tam_Vector_muestras);
        //printf("valor tmax - tmin %d\n",T_MAX-T_MIN);
        //printf("valor paso temporal %f\n", (double)(T_MAX-T_MIN) / Tam_Vector_muestras);
        double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;
        //printf("Valor paso temporal: %f\n",paso_temporal);
        clock_t inicio_cont = clock();
        CFT(Fourier_cont,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont = clock();
        double tiempo_cont = (double)(fin_cont-inicio_cont)/CLOCKS_PER_SEC;
        fprintf(salida_cont,"%d %lf\n",Tam_Vector_muestras,tiempo_cont);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont,"%lf %lf\n",creal(Fourier_cont[i]),cimag(Fourier_cont[i]));
        }*/
        //printf("Valor paso temporal: %f\n",paso_temporal);
        clock_t inicio_cont2 = clock();
        CFT_Simpson(Fourier_cont2,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont2 = clock();
        double tiempo_cont2 = (double)(fin_cont2-inicio_cont2)/CLOCKS_PER_SEC;
        fprintf(salida_cont2,"%d %lf\n",Tam_Vector_muestras,tiempo_cont2);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont2,"%lf %lf\n",creal(Fourier_cont2[i]),cimag(Fourier_cont2[i]));
        }*/
        //printf("Valor paso temporal: %f\n",paso_temporal);
        clock_t inicio_cont3 = clock();
        CFT_Trapecio(Fourier_cont3,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont3 = clock();
        double tiempo_cont3 = (double)(fin_cont3-inicio_cont3)/CLOCKS_PER_SEC;
        fprintf(salida_cont3,"%d %lf\n",Tam_Vector_muestras,tiempo_cont3);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont3,"%lf %lf\n",creal(Fourier_cont3[i]),cimag(Fourier_cont3[i]));
        }*/
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        free(muestras);
        free(Fourier);
        free(Fourier_cont);
        free(Fourier_cont2);
        free(Fourier_cont3);

    }

    fclose(entrada);
    fclose(salida);
    fclose(salida_cont);
    fclose(salida_cont2);
    fclose(salida_cont3);

    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier Discreto) en salidaDFT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo Simpson) en salidaCFT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo Trapezio) en salidaCFT\n");

    return 0;
    






}
