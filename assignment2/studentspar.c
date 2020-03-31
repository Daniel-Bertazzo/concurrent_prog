/*  Trabalho 1 - Programacao Concorrente - SSC0143
    Segundo semestre de 2019

    Daniel Penna Chaves Bertazzo - 10349561
    Gabriel Seiji Matsumoto      - 10295332
    Leonardo Miassi Netto        - 9326688

    PARA COMPILAR:
    mpicc studentspar.c -o studentspar -lm -fopenmp

    PARA RODAR:
    mpirun --oversubscribe -np N studentspar < input.in
    
    onde N e' o numero de processos (N >= 5) e input.in e' o arquivo
    que contem os dados a serem lidos.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include <omp.h>

// Struct que representa as regioes
typedef struct Regioes {
    int R, C, A; // Numero de regioes, cidades e alunos por cidade
    int *m; // Vetor que representa os dados da regiao
} Regioes;


/* ..:: Quicksort retirado de https://www.geeksforgeeks.org/quick-sort/ ::.. */
/********************************************************************************************/

void swap(int* a, int* b) { 
    int t = *a; 
    *a = *b; 
    *b = t; 
} 
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (int arr[], int low, int high) { 
    int pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j] < pivot) 
        { 
            i++;    // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 
  
/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
void quickSort(int arr[], int low, int high) { 
    if (low < high) 
    {   
        int pi = partition(arr, low, high); 
        
        #pragma omp task default(none) firstprivate(arr, low, pi)
        {
            quickSort(arr, low, pi - 1);
        }
        #pragma omp task default(none) firstprivate(arr, pi, high)
        {
            quickSort(arr, pi + 1, high); 
        }    
    } 
}
/********************************************************************************************/
/* ..:: TRABALHO ::.. */

// Le as dimensoes dos dados 
Regioes *le_dimensoes() {
    int R, A, C, seed;
    scanf("%d %d %d %d", &R, &C, &A, &seed);
    srand(seed);

    // Aloca as regioes
    Regioes *r = (Regioes *) malloc(sizeof(Regioes));
    
    r->R = R;
    r->C = C;
    r->A = A;
    r->m = NULL;
    
    return r;
}

// Monta matriz de dimensoes R*C x A (representada por um vetor)
// Essa unica matriz vai representar todas as regioes, com seus respectivos dados
void le_matriz(Regioes *r) {
    for (int i = 0; i < r->R*r->C; i++) {    
        for (int j = 0; j < r->A; j++) {
            r->m[i*r->A + j] = rand() % 101; // Valores de 0 a 100
        }
    }
}

/**********************************************************************************************************/
/* ..:: CALCULOS POR CIDADE ::.. */

// Calcula as maiores notas por cidade
void maiorCidade(Regioes *r, int *maioresCidade) {
    int i;

    // Vetor (r->m[i]) ja esta ordenado => maior nota na ultima posicao
    #pragma omp parallel for
    for (i = 0; i < r->R*r->C; i++) {
        maioresCidade[i] = r->m[(i*r->A) + (r->A-1)];
    }
}

// Calcula as menores notas por cidade
void menorCidade(Regioes *r, int *menoresCidade) {
    int i; 

    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    #pragma omp parallel for
    for (i = 0; i < r->R*r->C; i++) {
        menoresCidade[i] = r->m[i*r->A];
    }
}

// Calcula media aritmetica para cada cidade (entre os alunos)
void mediaAritmeticaCidade(Regioes *r, double *maCidade) {
    int i, j;

    #pragma omp parallel for private(i, j)
    for (i = 0; i < r->C*r->R; i++) {
        for (j = 0; j < r->A; j++) {
            maCidade[i] += r->m[i*r->A + j];
        }
        maCidade[i] = maCidade[i] / (double)r->A;
    }
}

// Calcula mediana para cada cidade (entre os alunos)
void medianaCidade(Regioes *r, double *medianasCidade){
    int i, j;
    int decisao, mid;

    // meio = num_colunas / 2
    mid = r->A/2;

    if((r->A%2) == 0) decisao = 0;
    else decisao = 1;

    #pragma omp parallel for
    for(i = 0; i < r->R * r->C; i++){
        if(decisao){
            medianasCidade[i] = r->m[i*r->A + mid];
        }
        else {
            medianasCidade[i] = (double)(r->m[i*r->A + (mid-1)] + r->m[(i*r->A) + (mid)]) / 2.0;
        }
    }
}

// Calcula o desvio padrao para cada cidade (entre os alunos)
void *desvioPadraoCidade(Regioes *r, double *dpCidade, double *maCidade) {
    int i, j;
    double n = r->A;

    #pragma omp parallel for private(i, j)
    for(i = 0; i < r->R * r->C; i++) {
        for(j = 0; j < r->A; j++) {
            dpCidade[i] += (r->m[i*r->A + j] - maCidade[i]) * (r->m[i*r->A + j] - maCidade[i]);
        }
        dpCidade[i] = sqrt(dpCidade[i]/(n-1.0));
    }
}

/**********************************************************************************************************/
/* ..:: CALCULOS POR REGIAO ::.. */

// Calcula as maiores notas por regiao
void maiorRegiao(Regioes *r, int *maioresRegiao) {
    // Vetor (r->m[i]) ja esta ordenado => maior nota na ultima posicao
    int tamRegiao = r->C * r->A;
    int i;

    #pragma omp parallel for
    for (i = 0; i < r->R; i++) {
        maioresRegiao[i] = r->m[(i*tamRegiao) + (tamRegiao-1)];
    }
}

// Calcula as menores notas por regiao
void menorRegiao(Regioes *r, int *menoresRegiao) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    int tamRegiao = r->C * r->A;
    int i;

    #pragma omp parallel for
    for (i = 0; i < r->R; i++) {
        menoresRegiao[i] = r->m[i*tamRegiao];
    }
}

// Calcula a media aritmetica para cada regiao (entre as cidades)
void mediaAritmeticaRegiao(Regioes *r, double *maRegiao) {
    int i, j;
    int tamRegiao = r->A * r->C;

    // Itera sobre as R regioes
    #pragma omp parallel for private(i, j)
    for (i = 0; i < r->R; i++) {
        // Itera sobre todos os dados de cada regiao
        for (j = 0; j < tamRegiao; j++) {
            maRegiao[i] += r->m[i*tamRegiao + j];
        }
        maRegiao[i] = maRegiao[i] / (double)tamRegiao;
    }
}

// Calcula mediana para cada Regiao 
void medianaRegiao(Regioes *r, double *medianasRegiao){
    int i, j;
    int decisao, mid;
    int tamRegiao = r->A * r->C;

    mid = (r->A * r->C)/2;

    if (((tamRegiao)%2) == 0) decisao = 0;
    else decisao = 1;

    #pragma omp parallel for 
    for(i = 0; i < r->R; i++){
        if(decisao) {
            medianasRegiao[i] = r->m[(i*tamRegiao) + (mid)];
        }
        else { 
            medianasRegiao[i] = (double)(r->m[i*tamRegiao + (mid-1)] + r->m[(i*tamRegiao) + (mid)]) / 2.0;
        }

    }
}

// Calcula o desvio padrao para cada regiao 
void *desvioPadraoRegiao(Regioes *r, double *dpRegiao, double *maRegiao) {
    int i,j;
    int tamRegiao = r->A * r->C;

    #pragma omp parallel for private(i, j)
    for(i = 0; i < r->R ; i++){
        for(j = 0; j < tamRegiao; j++) {
            dpRegiao[i] += (r->m[i*tamRegiao + j] - maRegiao[i]) * (r->m[i*tamRegiao + j] - maRegiao[i]);
        }
        dpRegiao[i] = sqrt(dpRegiao[i]/(double)(tamRegiao-1.0));
    }
}

/**********************************************************************************************************/
/* ..:: CALCULOS DO BRASIL ::.. */

// Calcula a maior nota do Brasil
int maiorBrasil(Regioes *r) {
    // Vetor (r->m) esta ordenado => maior na ultima posicao
    return r->m[(r->R * r->C * r->A) - 1];
}

// Calcula a menor nota do Brasil
int menorBrasil(Regioes *r) {
    // Vetor (r->m) esta ordenado => maior na primeira posicao
    return r->m[0];
}

// Calcula a media aritmetica para o Brasil
double mediaAritmeticaBrasil(Regioes *r) {
    double maBrasil = 0.0;
    int tamBrasil = r->R * r->C * r->A;
    int i;

    #pragma omp parallel for private(i) reduction(+:maBrasil)
    for (i = 0; i < tamBrasil; i++) {
        maBrasil += r->m[i];
    }

    return (maBrasil / (double)(tamBrasil));
}

// Calcula mediana para o Brasil
double medianaBrasil(Regioes *r) {
    int tamBrasil = r->R * r->C * r->A;
    int mid = tamBrasil / 2;

    if (((r->R * r->C * r->A) % 2) == 0) {
        return (double)(r->m[mid-1] + r->m[mid]) / 2.0;
    }
    else {
        return (double)r->m[mid];
    }

}

// Calcula o desvio padrao para  o Brasil
double desvioPadraoBrasil(Regioes *r, double maBrasil) {
    double soma = 0.0;
    double n = r->R*r->C*r->A;
    int i;

    #pragma omp parallel for private(i) shared(n) reduction(+:soma)
    for (i = 0; i < (int)n; i++){
        soma += (r->m[i] - maBrasil) * (r->m[i] - maBrasil);
    }
    return sqrt(soma/(n-1.0));
}


/* ..:: Exibição ::.. */
/********************************************************************************************/
void exibe(Regioes *r) {
    printf("R = %d, C = %d, A = %d\n", r->R, r->C, r->A);
    int lin = r->R * r->C;
    int col = r->A;

    for (int i = 0; i < lin; i++) {
        for (int j = 0; j < col; j++) {
            printf("%d ", r->m[i*r->A + j]);
        }
        printf("\n");
    }
}


/* ..:: MAIN ::.. */
/********************************************************************************************/

int main(int argc, char *argv[]) {
    int my_rank, num_procs, src, dest, msgtag;

    MPI_Status status;
    MPI_Request mpirequest_mr;
    Regioes *regioes = NULL;
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    // Se for o processo 0, le as dimensoes
    if (my_rank == 0) {
        regioes = le_dimensoes();
    }
    // Caso contrario, apenas aloca a struct
    else {
        regioes = (Regioes *) calloc(1, sizeof(Regioes));
        regioes->m = NULL;
    }

    // Manda as dimensoes para todos os outros processos
    src = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&regioes->R, 1, MPI_INT, src, MPI_COMM_WORLD);
    MPI_Bcast(&regioes->C, 1, MPI_INT, src, MPI_COMM_WORLD);
    MPI_Bcast(&regioes->A, 1, MPI_INT, src, MPI_COMM_WORLD);

    // printf("Existem %d processos. Eu sou o %d. R: %d, C: %d, A: %d\n", num_procs, my_rank, regioes->R, regioes->C, regioes->A);

    // Aloca o vetor para todos os processos
    regioes->m = (int *) calloc(regioes->R*regioes->C*regioes->A, sizeof(int));

    // Le os dados apenas no processo 0
    if (my_rank == 0) {
        le_matriz(regioes);
    }
    // Manda os dados lidos para os outros processos
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(regioes->m, regioes->R*regioes->C*regioes->A, MPI_INT, src, MPI_COMM_WORLD);

    int tamCidade = regioes->A;
    int tamRegiao = regioes->C * regioes->A;
    
    // Vetores que armazenam os dados para cada cidade
    int *maioresCidade     = (int *)    calloc(regioes->R * regioes->C, sizeof(int));
    int *menoresCidade     = (int *)    calloc(regioes->R * regioes->C, sizeof(int));
    double *maCidade       = (double *) calloc(regioes->R * regioes->C, sizeof(double)); // maCidade = medias aritmeticas das cidades
    double *medianasCidade = (double *) calloc(regioes->R * regioes->C, sizeof(double));
    double *dpCidade       = (double *) calloc(regioes->R * regioes->C, sizeof(double)); // dpCidade = desvios padroes das cidades

    // Vetores que armazenam os dados para cada regiao
    int *maioresRegiao     = (int *)    calloc(regioes->R, sizeof(int));
    int *menoresRegiao     = (int *)    calloc(regioes->R, sizeof(int));
    double *maRegiao       = (double *) calloc(regioes->R, sizeof(double));
    double *medianasRegiao = (double *) calloc(regioes->R, sizeof(double));
    double *dpRegiao       = (double *) calloc(regioes->R, sizeof(double));

    // Variaveis que armazenam os dados para o Brasil todo
    int majorBrasil = 0; 
    int minorBrasil = 0;
    double maBrasil = 0.0;
    double meBrasil = 0.0;
    double dpBrasil = 0.0;

    double t1 = MPI_Wtime();

    /* ..:: Calculos por cidade ::.. */
    if (my_rank == 1) {
        
        // Ordena as notas de cada cidade (ordena as linhas)
        #pragma omp parallel for shared(tamCidade)
        for (int i = 0; i < regioes->C*regioes->R; i++) {
            quickSort(regioes->m, i * tamCidade, (i*tamCidade) + (tamCidade-1));
        }

        // Chamando funcoes para as cidades
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                maiorCidade(regioes, maioresCidade);
            }
            #pragma omp section
            {
                menorCidade(regioes, menoresCidade);
            }
            #pragma omp section
            {
                mediaAritmeticaCidade(regioes, maCidade);
            }
            #pragma omp section
            {
                medianaCidade(regioes, medianasCidade);
            }
        }
        desvioPadraoCidade(regioes, dpCidade, maCidade);

        // Envia os resultados para o processo de impressao
        dest = 4;
        msgtag = 0;
        MPI_Send(maioresCidade,  regioes->R * regioes->C, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(menoresCidade,  regioes->R * regioes->C, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(maCidade,       regioes->R * regioes->C, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(medianasCidade, regioes->R * regioes->C, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(dpCidade,       regioes->R * regioes->C, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);

    }

    /* ..:: Calculos por regiao ::.. */
    else if (my_rank == 2) {
    
        // Ordena as notas das regioes (ordena os blocos CxA que representam as regioes)
        #pragma omp parallel for shared(tamRegiao)
        for (int i = 0; i < regioes->R; i++) {
            quickSort(regioes->m, i * tamRegiao, (i*tamRegiao) + (tamRegiao-1));
        }

        // Chamando funcoes para as regioes
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                maiorRegiao(regioes, maioresRegiao);    
            }
            #pragma omp section
            {
                menorRegiao(regioes, menoresRegiao);    
            }
            #pragma omp section
            {
                mediaAritmeticaRegiao(regioes, maRegiao);    
            }
            #pragma omp section
            {
                medianaRegiao(regioes, medianasRegiao);
            }
        }
        desvioPadraoRegiao(regioes, dpRegiao, maRegiao);

        // Envia os resultados para o processo de impressao
        dest = 4;
        msgtag = 0;
        MPI_Send(maioresRegiao,  regioes->R, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(menoresRegiao,  regioes->R, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(maRegiao,       regioes->R, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(medianasRegiao, regioes->R, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(dpRegiao,       regioes->R, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
    }

    /* ..:: Calculos para o Brasil ::.. */
    else if(my_rank == 3){  
        // Ordena todos os dados do pais
        int *a = regioes->m;
        int t = ( (regioes->R * regioes->C * regioes->A) - 1);
        #pragma omp parallel default (none) shared (a,t)
        {
            #pragma omp single nowait
            {
                quickSort(a, 0, t);
            }
        }

        // Chamando funcoes para o Brasil
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                majorBrasil = maiorBrasil(regioes);
            }
            #pragma omp section
            {
                minorBrasil = menorBrasil(regioes);
            }
            #pragma omp section
            {
                maBrasil = mediaAritmeticaBrasil(regioes);
            }
            #pragma omp section
            {
                meBrasil = medianaBrasil(regioes);
            }
        }
        dpBrasil = desvioPadraoBrasil(regioes,maBrasil);

        // Envia os resultados para o processo de impressao
        dest = 4;
        msgtag = 0;
        MPI_Send(&majorBrasil,  1, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(&minorBrasil,  1, MPI_INT,    dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(&maBrasil,     1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(&meBrasil,     1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
        MPI_Send(&dpBrasil,     1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
    }


    /* ..:: Imprimir os resultados ::.. */
    else if(my_rank == 4) {
    
        msgtag = 0;

        // Resultados das cidades
        src = 1;
        MPI_Recv(maioresCidade,  regioes->R * regioes->C, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(menoresCidade,  regioes->R * regioes->C, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(maCidade,       regioes->R * regioes->C, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(medianasCidade, regioes->R * regioes->C, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(dpCidade,       regioes->R * regioes->C, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);

        // Resultados das regioes
        src = 2;
        MPI_Recv(maioresRegiao,  regioes->R, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(menoresRegiao,  regioes->R, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(maRegiao,       regioes->R, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(medianasRegiao, regioes->R, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(dpRegiao,       regioes->R, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);

        // Resultados do Brasil
        src = 3;
        MPI_Recv(&majorBrasil,  1, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(&minorBrasil,  1, MPI_INT,    src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(&maBrasil,     1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(&meBrasil,     1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        MPI_Recv(&dpBrasil,     1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);

        double t2 = MPI_Wtime();

        double maiorMediaCidade = maCidade[0], maiorMediaRegiao = maRegiao[0];
        int qualCidade = 0, qualCidadeRegiao = 0, qualRegiao = 0;
        int n = regioes->A, cont = 0;


        for(int i = 0; i < regioes->R; i++){
            for(int j = 0; j < regioes->C; j++){
                printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2lf, media: %.2lf e DP: %.2lf\n",
                        i, j, menoresCidade[cont], maioresCidade[cont], medianasCidade[cont], maCidade[cont], dpCidade[cont]);

                if(maiorMediaCidade < maCidade[cont]){
                    maiorMediaCidade = maCidade[cont];
                    qualCidade = j;
                    qualCidadeRegiao = i;
                }
                cont++;
            }
        printf("\n");
        }
    
        //Resultado por Regiao
        
        for(int i = 0; i < regioes->R; i++){
            printf("Reg %d: menor: %d, maior: %d, mediana: %.2lf, media: %.2lf e DP: %.2lf\n",
                    i, menoresRegiao[i], maioresRegiao[i], medianasRegiao[i], maRegiao[i], dpRegiao[i]);    

            if(maiorMediaRegiao < maRegiao[i]){
                maiorMediaRegiao = maRegiao[i];
                qualRegiao = i;
            }
        }
        printf("\n");

        printf("Brasil: menor: %d, maior: %d, mediana: %.2lf, media: %.2lf e DP: %.2lf\n",
                minorBrasil, majorBrasil, meBrasil, maBrasil, dpBrasil);

        printf("\n");
        //Melhores Classificações
        printf("Melhor região: Região %d\n", qualRegiao);
        printf("Melhor cidade: Região %d, Cidade %d\n", qualCidadeRegiao, qualCidade);

        printf("Tempo de resposta sem considerar E/S, em segundos: %lfs\n", t2-t1);
    }


    MPI_Finalize();

    /* ..:: Liberacao de memoria ::.. */
    free(regioes->m);
    free(regioes);
    free(maioresCidade);
    free(menoresCidade);
    free(maCidade);
    free(medianasCidade);
    free(dpCidade);
    free(maioresRegiao);
    free(menoresRegiao);
    free(maRegiao);
    free(medianasRegiao);
    free(dpRegiao);

    return 0;
}