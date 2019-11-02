/*  Trabalho 1 - Programacao Concorrente - SSC0143
    Segundo semestre de 2019

    Daniel Penna Chaves Bertazzo - 10349561
    Gabriel Seiji Matsumoto      - 10295332
    Leonardo Miassi Netto        - 9326688

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Struct que representa as regioes
typedef struct Regioes {
    int R, C, A; // Numero de regioes, cidades e alunos por cidade
    int *m; // Vetor que representa os dados da regiao
} Regioes;


/* ..:: QUICKSORT RETIRADO DO GEEKS FOR GEEKS ::.. */
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
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
}

/* ..:: TRABALHO ::.. */

Regioes *le_entrada() {
    // Le a entrada do arquivo
    int R, A, C, seed;
    int i, j;
    scanf("%d %d %d %d", &R, &C, &A, &seed);

    // Aloca as regioes e atribui os valores
    Regioes *aux = (Regioes *) malloc(sizeof(Regioes));
    srand(seed);
    aux->R = R;
    aux->C = C;
    aux->A = A;

    // Monta matriz de dimensoes R*C x A (representada por um vetor)
    // Essa unica matriz vai representar todas as regioes, com seus respectivos dados
    aux->m = (int *) malloc(R*C*A * sizeof(int *));
    for(i = 0; i < R*C; i++) {    
        for(j = 0; j < A; j++) {
            aux->m[i*A + j] = rand() % 100;
        }
    }
    
    return aux;
}

/**********************************************************************************************************/
/* ..:: CALCULOS POR CIDADE ::.. */

// Calcula as maiores notas por cidade
void maiorCidade(Regioes *r, int *maioresCidade) {
    // Vetor (r->m[i]) ja esta ordenado => maior nota na ultima posicao
    for (int i = 0; i < r->R*r->C; i++) {
        // maioresCidade[i] = r->m[i][r->A-1];
        maioresCidade[i] = r->m[(i*r->A) + (r->A-1)];
    }
}


// Calcula as menores notas por cidade
void menorCidade(Regioes *r, int *menoresCidade) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    for (int i = 0; i < r->R*r->C; i++) {
        // menoresCidade[i] = r->m[i][0];
        menoresCidade[i] = r->m[i*r->A];
    }
}


// Calcula media aritmetica para cada cidade (entre os alunos)
void MediaAritmeticaCidade(Regioes *r, double *maCidade) {
    int i,j;

    for (i = 0; i < r->C*r->R; i++) {
        for (j = 0; j < r->A; j++) {
            // maCidade[i] += r->m[i][j];
            maCidade[i] += r->m[i*r->A + j];
        }
        maCidade[i] = maCidade[i] / (double)r->A;
    }
}

// Calcula mediana para cada cidade (entre os alunos)
void MedianaCidade(Regioes *r, double *medianasCidade){
    int i,j;
    int decisao, mid;

    mid = r->A/2;

    if((r->A%2) == 0) decisao = 0;
    else decisao = 1;

    for(i = 0; i < r->R * r->C; i++){
        if(decisao)
            // medianasCidade[i] = r->m[i][mid+1];
            medianasCidade[i] = r->m[(i*r->A) + (mid+1)];
        else 
            // medianasCidade[i] = (r->m[i][mid] + r->m[i][mid+1]) / 2.0;
            medianasCidade[i] = (r->m[i*r->A + mid] + r->m[(i+r->A) + (mid+1)]) / 2.0;

    }
}

// Calcula o desvio padrao para cada cidade (entre os alunos)
void *desvioPadraoCidade(Regioes *r, double *dpCidades, double *maCidade) {
    int i, j;

    for(i = 0; i < r->R * r->C; i++) {
        for(j = 0; j < r->A; j++) {
            // dpCidades[i] += (r->m[i][j] - maCidade[i]) * (r->m[i][j] - maCidade[i]);
            dpCidades[i] += (r->m[i*r->A + j] - maCidade[i]) * (r->m[i*r->A + j] - maCidade[i]);
        }
        dpCidades[i] = sqrt(dpCidades[i]/(double)(r->R * r->C));
    }
}

/**********************************************************************************************************/
/* ..:: CALCULOS POR REGIAO ::.. */

// Calcula as maiores notas por regiao
void maiorRegiao(Regioes *r, int *maioresRegiao) {
    // Vetor (r->m[i]) ja esta ordenado => maior nota na ultima posicao
    int tamRegiao = r->C * r->A;
    for (int i = 0; i < r->R; i++) {
        // maioresRegiao[i] = r->m[i][r->A-1];
        maioresRegiao[i] = r->m[(i*tamRegiao) + (tamRegiao-1)];
    }
}


// Calcula as menores notas por regiao
void menorRegiao(Regioes *r, int *menoresRegiao) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    int tamRegiao = r->C * r->A;
    for (int i = 0; i < r->R; i++) {
        // menoresRegiao[i] = r->m[i][0];
        menoresRegiao[i] = r->m[i*tamRegiao];
    }
}

// Calcula a media aritmetica para cada regiao (entre as cidades)
void MediaAritmeticaRegiao(Regioes *r, double *maRegiao) {
    int i, j;
    int tamRegiao = r->A * r->C;

    // Itera sobre as R regioes
    for (i = 0; i < r->R; i++) {
        // Itera sobre todos os dados de cada regiao
        for (j = i * tamRegiao; j < (i*tamRegiao) + (tamRegiao-1); j++) {
            maRegiao[i] += r->m[i*r->A + j];
        }
        maRegiao[i] = maRegiao[i] / (double)tamRegiao;
    }
}

// Calcula mediana para cada Regiao 
void MedianaRegiao(Regioes *r, double *medianasRegiao){
    int i, j;
    int decisao, mid;

    mid = (r->A * r->C)/2;

    if (((r->A * r->C)%2) == 0) decisao = 0;
    else decisao = 1;

    for(i = 0; i < r->R; i++){
        if(decisao)
            // medianasRegiao[i] = r->m[i][mid+1];
            medianasRegiao[i] = r->m[(i*r->A) + (mid+1)];
        else 
            // medianasRegiao[i] = (r->m[i][mid] + r->m[i][mid+1]) / 2.0;
            medianasRegiao[i] = (r->m[i*r->A + mid] + r->m[(i+r->A) + (mid+1)]) / 2.0;

    }
}

// Calcula o desvio padrao para cada regiao 
void *desvioPadraoRegiao(Regioes *r, double *dpRegiao, int *maRegiao) {
    int i,j;
    int tamRegiao = r->A * r->C;

    for(i = 0; i < r->R ; i++){
        for(j = i * tamRegiao; j < (i*tamRegiao) + (tamRegiao-1); j++) {
            // dpRegiao[i] += (r->m[i][j] - medias[i]) * (r->m[i][j] - medias[i]);
            dpRegiao[i] += (r->m[i*r->A + j] - maRegiao[i]) * (r->m[i*r->A + j] - maRegiao[i]);
        }
        dpRegiao[i] = sqrt(dpRegiao[i]/(double)(tamRegiao));
    }
}

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

/* ..:: Liberacao de memoria ::.. */
/********************************************************************************************/
void libera_memoria(Regioes *r) {
    free(r->m);
    free(r);
}


/* ..:: MAIN ::.. */
/********************************************************************************************/

int main(int argc, char const *argv[]) {
    int i, j;
    Regioes *regioes = le_entrada();

    /* ..:: Alocacao dos dados necessarios ::.. */
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
    

    // Ordena as notas de cada cidade (ordena as linhas)
    int tamCidade = regioes->A;
    for (i = 0; i < regioes->C*regioes->R; i++) {
        quickSort(regioes->m, i * tamCidade, (i*tamCidade) + (tamCidade-1));
    }

    /* ..:: Chamando funcoes para as cidades ::.. */



    // Ordena as notas das regioes (ordena os blocos CxA que representam as regioes)
    int tamRegiao = regioes->C * regioes->A;
    for (i = 0; i < regioes->R; i++) {
        quickSort(regioes->m, i * tamRegiao, (i*tamRegiao) + (tamRegiao-1));
    }




    /* ..:: Liberacao de memoria ::.. */
    libera_memoria(regioes);
    
    return 0;
}