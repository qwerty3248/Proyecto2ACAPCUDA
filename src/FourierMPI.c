#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>

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
/*void DFT(double complex *Fourier,const double *muestras, const int n){
    for (int i=0;i<n;i++){
        Fourier[i] = 0.0 + 0.0*I;
        for (int j=0;j<n;j++){
            double angle = -2.0*PI*i*j/n;
            Fourier[i] += muestras[j]*cexp(I*angle);
        }
    }    
}*/


//Fourier continuo 🤯🤯
/*void CFT(double complex *Fourier, const double *muestras, const int n, const double paso_temporal){
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
}*/

//Integral con simpson (2)

/*void CFT_Simpson(double complex *Fourier, const double *muestras, const int n, const double paso_temporal) {
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
}*/


//Integral con trapecio (3)
/*void CFT_Trapecio(double complex *Fourier, const double *muestras, const int n, const double paso_temporal) {
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
}*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Aqui las funciones de nuevo pero para paralelizar 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DFT(double complex *Fourier, const double *muestras, const int n, const int rank, const int size) {
    int local = n / size;
    if (rank < n % size) {
        local++;
    }
    int comienzo = rank * local;
    int final = (rank == size - 1) ? n : comienzo + local;

    // Reserva memoria para la copia local de Fourier
    double complex *copia = (double complex *)malloc(local * sizeof(double complex));

    if (!copia) {
        printf("Error: No se pudo asignar memoria para copia\n");
        exit(1);
    }

    // Calcula la DFT para la copia local
    for (int i = comienzo; i < final; i++) {
        copia[i - comienzo] = 0.0 + 0.0 * I;
        for (int j = 0; j < n; j++) {
            double angle = -2.0 * PI * i * j / n;
            copia[i - comienzo] += muestras[j] * cexp(I * angle);
        }
    }

    // Comunica la copia local a la memoria compartida
    MPI_Gather(copia, local, MPI_DOUBLE_COMPLEX, Fourier, local, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    
    free(copia);
}


//Fourier continuo 🤯🤯 Riemann
void CFT(double complex *Fourier, const double *muestras, const int n, const double paso_temporal, const int rank, const int size) {
    int local = n / size;
    if (rank < n % size) {
        local++;
    }
    int comienzo = rank * local;
    int final = (rank == size - 1) ? n : comienzo + local;

    
    double complex *copia = (double complex *)malloc(sizeof(double complex) * local);

    if (!copia) {
        printf("Error: No se pudo asignar memoria para copia\n");
        exit(1);
    }
    
    for (int i = comienzo; i < final; i++) {
        copia[i - comienzo] = 0.0 + 0.0 * I;
        double omega = 2.0 * PI * i / (T_MAX - T_MIN);
        for (int j = 0; j < n; j++) {
            double t = T_MIN + j * paso_temporal;
            if (t < T_MIN) {
                copia[i - comienzo] += muestras[0] * cexp(-I * omega * t) * paso_temporal;
            } else if (t >= T_MAX) {
                copia[i - comienzo] += muestras[n - 1] * cexp(-I * omega * t) * paso_temporal;
            } else {
                copia[i - comienzo] += muestras[j] * cexp(-I * omega * t) * paso_temporal;
            }
        }
    }

    
    MPI_Gather(copia, local, MPI_DOUBLE_COMPLEX, Fourier, local, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    
    free(copia);
}

void CFT_Trapecio(double complex *Fourier, const double *muestras, const int n, const double paso_temporal,const int rank, const int size) {
    int local = n / size;
    if (rank < n % size) {
        local++;
    }
    int comienzo = rank * local;
    int final = (rank == size - 1) ? n : comienzo + local;    
    double complex *copia = (double complex *)malloc(sizeof(double complex) * local);
    if (!copia) {
        printf("Error: No se pudo asignar memoria para copia\n");
        exit(1);
    }
    for (int i = comienzo; i < final; i++) {
        copia[i - comienzo] = 0.0 + 0.0 * I;
        double omega = 2.0 * PI * i / (T_MAX - T_MIN);
        for (double j = T_MIN; j < T_MAX - paso_temporal; j += paso_temporal) {
            int indice1 = (int)((j - T_MIN) / paso_temporal);
            int indice2 = indice1 + 1;
            if (indice2 >= n) indice2 = n - 1;
            double promedio = (muestras[indice1] + muestras[indice2]) / 2.0;
            copia[i - comienzo] += promedio * cexp(-I * omega * j) * paso_temporal;
        }
    }

    MPI_Gather(copia, local, MPI_DOUBLE_COMPLEX, Fourier, local, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    free(copia);

}

void CFT_Simpson(double complex *Fourier, const double *muestras, const int n, const double paso_temporal,const int rank, const int size) {
    int local = n / size;
    if (rank < n % size) {
        local++;
    }
    int comienzo = rank * local;
    int final = (rank == size - 1) ? n : comienzo + local;

    
    double complex *copia = (double complex *)malloc(sizeof(double complex) * local);

    if (!copia) {
        printf("Error: No se pudo asignar memoria para copia\n");
        exit(1);
    }

    for (int i = comienzo; i < final; i++) {
        copia[i - comienzo] = 0.0 + 0.0 * I;
        double omega = 2.0 * PI * i / (T_MAX - T_MIN);
        for (double j = T_MIN; j < T_MAX - paso_temporal; j += paso_temporal) {
            int indice1 = (int)((j - T_MIN) / paso_temporal);
            int indice2 = indice1 + 1;

            if (indice2 >= n) indice2 = n - 1;

            double x_medio = j + paso_temporal / 2.0;
            int indice_medio = (int)((x_medio - T_MIN) / paso_temporal);
            if (indice_medio >= n) indice_medio = n - 1;

            double simpson = (muestras[indice1] + 4.0 * muestras[indice_medio] + muestras[indice2]) / 6.0;

            copia[i - comienzo] += simpson * cexp(-I * omega * j) * paso_temporal;
        }
    }
    //Comunica la copia local a la memoria compartida
    MPI_Gather(copia, local, MPI_DOUBLE_COMPLEX, Fourier, local, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    
    free(copia);
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
const char * SALIDA = "txt/SecuencialDFTMPI.txt"; //salida secuencial (cualquier caso) 😇😇
const char * SALIDA2 = "txt/SecuencialDFT2MPI.txt"; //salida secuencial para el caso de 20.000 muestras 😇😇
const char * SALIDA_CONTINUO = "txt/ContinuoDFTMPI.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "txt/ContinuoDFTMPI2.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO3 = "txt/ContinuoDFTMPI3.txt"; //salida continuo 🤯🤯
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (){
    
    MPI_Init(NULL,NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    FILE * salida;
    FILE * salida_cont;
    FILE * salida_cont2;
    FILE * salida_cont3;
    if (rank == 0){
        //printf("Archivo de salida: %s\n",new_salida);
        salida = fopen(SALIDA,"w");
        salida_cont = fopen(SALIDA_CONTINUO,"w");
        salida_cont2 = fopen(SALIDA_CONTINUO2,"w");
        salida_cont3 = fopen(SALIDA_CONTINUO3,"w");
    
        
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



    }//Solo quiero que reserve memoria en el 0, este se encargara de escribir, este tendra el fourier
    //Todos abren el archivo de entrada
    FILE * entrada = fopen (DOSCIENTOSCINCUENTAMIL,"r");
    if (!entrada){
        printf("Error: No se pudo abrir el archivo de entrada\n");
        exit(1);
    } 

    //printf("Soy el procesador %d y voy a leer el archivo de entrada\n",rank);
    
    int Tam_Vector_muestras;
    while (fscanf(entrada,"%d",&Tam_Vector_muestras) == 1){
        
        double *muestras = malloc(Tam_Vector_muestras*sizeof(double));
        double complex *Fourier = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont2 = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont3 = malloc(Tam_Vector_muestras*sizeof(double complex));
        
        
        //printf("Soy el procesador %d y voy a asignar memoria\n",rank);
        
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
        
        
        //Parte  del DFT 
        MPI_Barrier(MPI_COMM_WORLD);
        
        double tiempo_inicial = MPI_Wtime();
        DFT(Fourier, muestras,Tam_Vector_muestras,rank,size);
        double tiempo_final = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            fprintf(salida, "%d %lf\n", Tam_Vector_muestras, tiempo_final - tiempo_inicial);
        }
        /*if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida,"%lf %lf\n",creal(Fourier[i]),cimag(Fourier[i]));
            }
        }*/
        

        
        //double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;

        //Parte del CFT
        MPI_Barrier(MPI_COMM_WORLD);
        double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;
        tiempo_inicial = MPI_Wtime();
        CFT(Fourier_cont,muestras,Tam_Vector_muestras,paso_temporal,rank,size);
        tiempo_final = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            fprintf(salida_cont, "%d %lf\n", Tam_Vector_muestras, tiempo_final - tiempo_inicial);
        }
        /*if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida_cont,"%lf %lf\n",creal(Fourier_cont[i]),cimag(Fourier_cont[i]));
            }
        }*/
        //printf("andamos caminando\n");
        MPI_Barrier(MPI_COMM_WORLD);
        tiempo_inicial = MPI_Wtime();
        CFT_Trapecio(Fourier_cont3,muestras,Tam_Vector_muestras,paso_temporal,rank,size);
        tiempo_final = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            fprintf(salida_cont3, "%d %lf\n", Tam_Vector_muestras, tiempo_final - tiempo_inicial);
        }
        
        /*if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida_cont3,"%lf %lf\n",creal(Fourier_cont3[i]),cimag(Fourier_cont3[i]));
            }
        }*/

        //printf("andamos caminando2\n");
        MPI_Barrier(MPI_COMM_WORLD);
        tiempo_inicial = MPI_Wtime();
        CFT_Simpson(Fourier_cont2,muestras,Tam_Vector_muestras,paso_temporal,rank,size);
        tiempo_final = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            fprintf(salida_cont2, "%d %lf\n", Tam_Vector_muestras, tiempo_final - tiempo_inicial);
        }
        /*if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida_cont2,"%lf %lf\n",creal(Fourier_cont2[i]),cimag(Fourier_cont2[i]));
            }
        }*/
        
        
        

        //Todos los procesos liberan memoria
        free(muestras);
        free(Fourier);
        free(Fourier_cont);
        free(Fourier_cont2);
        free(Fourier_cont3);

    }

    //Solo lo cierra el proceso 0 porque es el que escribe 
    if (rank == 0){
        fclose(salida);
        fclose(salida_cont);
        fclose(salida_cont2);
        fclose(salida_cont3);
    }
    //Todos los procesos cierran el archivo de entrada porque todos los abren 
    fclose(entrada);
    //Para que todos los procesos terminen
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier Discreto) en salidaDFTMPI\n");
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFTMPI\n");
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFTMPI2\n");
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFTMPI3\n");
    }
    MPI_Finalize();
    return 0;
    
}