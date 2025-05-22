#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuComplex.h>


#define PI 3.14159265358979323846
//Las muestras están entre -1 y 1
#define T_MAX 1
#define T_MIN -1
#define BLOCK_SIZE 256
#define NUM_BLOCKS 10

//Si no hay exponencial en cuda la creo yo (Informacion cogida de la misma pagina de CUDA: https://forums.developer.nvidia.com/t/additional-cucomplex-functions-cucnorm-cucsqrt-cucexp-and-some-complex-double-functions/36892)
__host__ __device__ static __inline__ cuDoubleComplex cuCexp(cuDoubleComplex x)
{
	double factor = exp(x.x);
	return make_cuDoubleComplex(factor * cos(x.y), factor * sin(x.y));
}

//Esta si que la he puesto yo como me gusta, pero tambine viene en la pagina como un CuMul normal para los complejos dobles

__host__ __device__ static __inline__ cuDoubleComplex cuCmulReal(cuDoubleComplex a, double r) {
    return make_cuDoubleComplex(a.x * r, a.y * r);
}


//Aqui mis funciones 
__global__ void DFT(cuDoubleComplex *Fourier, const double *muestras, const int TAM_VECTOR_MUESTRAS){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < TAM_VECTOR_MUESTRAS){
        Fourier[i] = make_cuDoubleComplex(0.0,0.0);
        cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
        for (int j=0;j<TAM_VECTOR_MUESTRAS;j++){
            double angle = -2.0*PI*i*j/TAM_VECTOR_MUESTRAS;
            cuDoubleComplex aux = make_cuDoubleComplex(muestras[j]*cos(angle),muestras[j]*sin(angle));
            sum = cuCadd(sum,aux);
        }
        Fourier[i] = sum;
    }    
}

__global__ void CFT(cuDoubleComplex *Fourier, const double *muestras, const int TAM_VECTOR_MUESTRAS, const double paso_temporal){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < TAM_VECTOR_MUESTRAS){
        Fourier[i] = make_cuDoubleComplex(0.0,0.0);
        cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega
        //Aqui ya entra en juego tanto el intervalo del tiempo como los valores que dan la funcion, es decir el tiempo esta entre -1 y 1
        //y los valores de la funcion están en mis muestras
        //Vamos desde el minimo hasta el maximo pero con nuestro paso temporal para tomar fourier lo más preciso posible
        for (double j = T_MIN;j<T_MAX; j= j+paso_temporal){
            //printf("Estoy en el segundo FOR");
            int indice = (int)((j - T_MIN)/paso_temporal);
            if (indice <= 0 ){
                // Si es menor o igual a 0, suponemos que coge el primer elemento
                cuDoubleComplex expo = cuCexp(make_cuDoubleComplex(0.0, -omega * j));  
                cuDoubleComplex temp = make_cuDoubleComplex(muestras[0], 0.0);  
                cuDoubleComplex prod = cuCmul(temp, expo);  
                sum = cuCadd(sum, cuCmulReal(prod, paso_temporal));  
            }
            else if (indice >= TAM_VECTOR_MUESTRAS - 1) { 
                // Si es mayor o igual al número de elementos, cogemos el último
                cuDoubleComplex expo = cuCexp(make_cuDoubleComplex(0.0, -omega * j));  
                cuDoubleComplex temp = make_cuDoubleComplex(muestras[TAM_VECTOR_MUESTRAS - 1], 0.0);  
                cuDoubleComplex prod = cuCmul(temp, expo);  
                sum = cuCadd(sum, cuCmulReal(prod, paso_temporal));  
            } else {
                // Si no es válido y cogemos el valor calculado
                cuDoubleComplex expo = cuCexp(make_cuDoubleComplex(0.0, -omega * j));  // Usamos la función cuCexp para la exponencial
                cuDoubleComplex temp = make_cuDoubleComplex(muestras[indice], 0.0);  // Tomamos la muestra correspondiente
                cuDoubleComplex prod = cuCmul(temp, expo);  // Multiplicamos la muestra por la exponencial
                sum = cuCadd(sum, cuCmulReal(prod, paso_temporal));  // Acumulamos el resultado, aplicando el paso temporal
            }
        }
        Fourier[i] = sum;
    }    
}

__global__ void CFT_Simpson(cuDoubleComplex *Fourier, const double *muestras, const int TAM_VECTOR_MUESTRAS, const double paso_temporal){

        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < TAM_VECTOR_MUESTRAS){
            Fourier[i] = make_cuDoubleComplex(0.0,0.0);
            cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
            double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega

            for (double j = T_MIN;j<T_MAX - paso_temporal; j += paso_temporal){
                int indice1 = (int)((j - T_MIN)/paso_temporal);
                int indice2 = indice1 + 1;

                if (indice2 >= TAM_VECTOR_MUESTRAS) indice2 = TAM_VECTOR_MUESTRAS - 1;
                
                double x_medio = j + paso_temporal / 2.0;
                int indice_medio = (int)((x_medio - T_MIN) / paso_temporal);

                if (indice_medio >= TAM_VECTOR_MUESTRAS) indice_medio = TAM_VECTOR_MUESTRAS - 1;

                double simpson = (muestras[indice1] + 4.0*muestras[indice_medio] + muestras[indice2])/6.0;

                //Aqui la parte nueva además de distribuirlo para cada 
                cuDoubleComplex prod = cuCexp(make_cuDoubleComplex(0.0, -omega * j));
                cuDoubleComplex temp = make_cuDoubleComplex(simpson, 0.0);
                cuDoubleComplex prod2 = cuCmul(temp, prod);
                sum = cuCadd(sum,cuCmulReal(prod2,paso_temporal));


            }

            Fourier[i] = sum;

        }


}

__global__ void CFT_Trapecio(cuDoubleComplex *Fourier, const double *muestras, const int TAM_VECTOR_MUESTRAS, const double paso_temporal){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < TAM_VECTOR_MUESTRAS){
        Fourier[i] = make_cuDoubleComplex(0.0,0.0);
        cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega

        for(double j = T_MIN; j < T_MAX - paso_temporal; j += paso_temporal){
            int indice1 = (int)((j - T_MIN) / paso_temporal);
            int indice2 = indice1 + 1;

            if (indice2 >= TAM_VECTOR_MUESTRAS) indice2 = TAM_VECTOR_MUESTRAS - 1;

            double promedio = (muestras[indice1] + muestras[indice2])/2.0;

            //Aqui la parte nueva además de distribuirlo para cada
            cuDoubleComplex prod = cuCexp(make_cuDoubleComplex(0.0, -omega * j));
            cuDoubleComplex temp = make_cuDoubleComplex(promedio, 0.0);
            cuDoubleComplex prod2 = cuCmul(temp, prod);
            sum = cuCadd(sum,cuCmulReal(prod2,paso_temporal));

        }

        Fourier[i] = sum;

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
const char * SALIDA = "Resultados/SecuencialCUDADFT.txt"; //salida secuencial (cualquier caso) 😇😇
const char * SALIDA2 = "txt/SecuencialDFT2.txt"; //salida secuencial para el caso de 20.000 muestras 😇😇
const char * SALIDA_CONTINUO = "Resultados/SecuencialCUDACFTNormal.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "Resultados/SecuencialCUDACFTSimpson.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO3 = "Resultados/SecuencialCUDACFTTrapecio.txt"; //salida continuo 🤯🤯
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
    double *muestras;
    cuDoubleComplex *Fourier;
    cuDoubleComplex *Fourier_cont;
    cuDoubleComplex *Fourier_cont2;
    cuDoubleComplex *Fourier_cont3;
    while (fscanf(entrada,"%d",&Tam_Vector_muestras) == 1){
        
        cudaMallocManaged(&muestras,Tam_Vector_muestras*sizeof(double));//double *muestras = malloc(Tam_Vector_muestras*sizeof(double));
        cudaMallocManaged(&Fourier ,Tam_Vector_muestras*sizeof(cuDoubleComplex) );//double complex *Fourier = malloc(Tam_Vector_muestras*sizeof(double complex));
        cudaMallocManaged(&Fourier_cont ,Tam_Vector_muestras*sizeof(cuDoubleComplex) );//double complex *Fourier_cont = malloc(Tam_Vector_muestras*sizeof(double complex));
        cudaMallocManaged(&Fourier_cont2 ,Tam_Vector_muestras*sizeof(cuDoubleComplex) );//double complex *Fourier_cont2 = malloc(Tam_Vector_muestras*sizeof(double complex));
        cudaMallocManaged(&Fourier_cont3 ,Tam_Vector_muestras*sizeof(cuDoubleComplex) );//double complex *Fourier_cont3 = malloc(Tam_Vector_muestras*sizeof(double complex));
        
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

        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        cudaEventRecord(start, 0);
        DFT<<<NUM_BLOCKS,BLOCK_SIZE>>>(Fourier,muestras,Tam_Vector_muestras);
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        
        cudaEventElapsedTime(&milliseconds, start, stop);
        
        fprintf(salida,"%d %lf\n",Tam_Vector_muestras,milliseconds);
        //Para mostrar el vector, esta dentro del archivo 😴😴
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida, "%lf %lf\n", cuCreal(Fourier[i]), cuCimag(Fourier[i]));
        }*/

        cudaEventDestroy(start);
        cudaEventDestroy(stop);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de CFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //printf("valor tammuestras %d\n",Tam_Vector_muestras);
        //printf("valor tmax - tmin %d\n",T_MAX-T_MIN);
        //printf("valor paso temporal %f\n", (double)(T_MAX-T_MIN) / Tam_Vector_muestras);
        
        
        double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;
        cudaEvent_t start1, stop1;
        milliseconds = 0;
        cudaEventCreate(&start1);
        cudaEventCreate(&stop1);


        cudaEventRecord(start1, 0);
        CFT<<<NUM_BLOCKS,BLOCK_SIZE>>>(Fourier_cont,muestras,Tam_Vector_muestras,paso_temporal);
        cudaEventRecord(stop1, 0);
        cudaEventSynchronize(stop1);

        cudaEventElapsedTime(&milliseconds, start1, stop1);
        
        //printf("Valor paso temporal: %f\n",paso_temporal);
        //clock_t inicio_cont = clock();
        //CFT(Fourier_cont,muestras,Tam_Vector_muestras,paso_temporal);
        //clock_t fin_cont = clock();
        //double tiempo_cont = (double)(fin_cont-inicio_cont)/CLOCKS_PER_SEC;
        fprintf(salida_cont,"%d %lf\n",Tam_Vector_muestras,milliseconds);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont, "%lf %lf\n", cuCreal(Fourier_cont[i]), cuCimag(Fourier_cont[i]));
        }*/

        //printf("Valor paso temporal: %f\n",paso_temporal);
        /*clock_t inicio_cont2 = clock();
        CFT_Simpson(Fourier_cont2,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont2 = clock();
        double tiempo_cont2 = (double)(fin_cont2-inicio_cont2)/CLOCKS_PER_SEC;
        fprintf(salida_cont2,"%d %lf\n",Tam_Vector_muestras,tiempo_cont2);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont2,"%lf %lf\n",creal(Fourier_cont2[i]),cimag(Fourier_cont2[i]));
        }*/

        cudaEvent_t start2, stop2;
        float milliseconds2 = 0;
        cudaEventCreate(&start2);
        cudaEventCreate(&stop2);

        cudaEventRecord(start2, 0);
        CFT_Simpson<<<NUM_BLOCKS,BLOCK_SIZE>>>(Fourier_cont2,muestras,Tam_Vector_muestras,paso_temporal);
        cudaEventRecord(stop2, 0);
        cudaEventSynchronize(stop2);

        cudaEventElapsedTime(&milliseconds2, start2, stop2);
        
        fprintf(salida_cont2,"%d %lf\n",Tam_Vector_muestras,milliseconds2);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont2, "%lf %lf\n", cuCreal(Fourier_cont2[i]), cuCimag(Fourier_cont2[i]));
        }*/



        
        cudaEvent_t start3, stop3;
        float milliseconds3 = 0;
        cudaEventCreate(&start3);
        cudaEventCreate(&stop3);

        cudaEventRecord(start3, 0);
        CFT_Trapecio<<<NUM_BLOCKS,BLOCK_SIZE>>>(Fourier_cont3,muestras,Tam_Vector_muestras,paso_temporal);
        cudaEventRecord(stop3, 0);
        cudaEventSynchronize(stop3);

        cudaEventElapsedTime(&milliseconds3, start3, stop3);
        
        fprintf(salida_cont3,"%d %lf\n",Tam_Vector_muestras,milliseconds3);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont3, "%lf %lf\n", cuCreal(Fourier_cont3[i]), cuCimag(Fourier_cont3[i]));
        }*/

        


        //printf("Valor paso temporal: %f\n",paso_temporal);
        /*clock_t inicio_cont3 = clock();
        CFT_Trapecio(Fourier_cont3,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont3 = clock();
        double tiempo_cont3 = (double)(fin_cont3-inicio_cont3)/CLOCKS_PER_SEC;
        fprintf(salida_cont3,"%d %lf\n",Tam_Vector_muestras,tiempo_cont3);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont3,"%lf %lf\n",creal(Fourier_cont3[i]),cimag(Fourier_cont3[i]));
        }*/
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        cudaFree(muestras);
        cudaFree(Fourier);
        cudaFree(Fourier_cont);
        cudaFree(Fourier_cont2);
        cudaFree(Fourier_cont3);

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
