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
    int **m; // Matriz que representa os dados da regiao
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
/********************************************************************************************/

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

    // Monta matriz de dimensoes R*C x A
    // Essa unica matriz vai representar todas as regioes, com seus respectivos dados
    aux->m = (int **) malloc(R*C * sizeof(int *));
    for(i = 0; i < R*C; i++){    
        aux->m[i] = (int *) malloc(A * sizeof(int));
        for(j = 0; j < A; j++){
            aux->m[i][j] = rand() % 100;
        }
    }
    
    return aux;
}

// Calcula as maiores notas por cidade
void maiorCidade(Regioes *r, int *maiores) {
    // Vetor (r->m[i]) ja esta ordenado => maior nota na ultima posicao
    for (int i = 0; i < r->R*r->C; i++) {
        maiores[i] = r->m[i][A-1];
    }
}

// Calcula as menores notas por cidade
void menorCidade(Regioes *r, int *menores) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    for (int i = 0; i < r->R*r->C; i++) {
        menores[i] = r->m[i][0];
    }
}

// Calcula media aritmetica para cada cidade (entre os alunos)
void MediaAritmeticaCidade(Regioes *reg, double *maCidade) {
    int i,j;

    for (i = 0; i < reg->C*reg->R; i++) {
        for (j = 0; j < reg->A; j++) {
            maCidade[i] += reg->m[i][j];
        }
        maCidade[i] = maCidade[i] / (double)reg->A;
    }
}

// Calcula mediana para cada cidade (entre os alunos)
double *MedianaCidade(Regioes *reg, double *medianasCidade){
    int i,j;
    int Decisao,mid;

    mid = reg->A/2;

    if((reg->A%2) == 0) Decisao = 0;
    else Decisao = 1;

    for(i = 0; i < reg->R * reg->C; i++){
        if(Decisao)
            Medianas[i] = reg->m[i][mid+1];
        else 
            Medianas[i] = (reg->m[i][mid] + reg->m[i][mid+1]) / 2.0;
    }

    return vet;
}

// Calcula o desvio padrao para cada cidade (entre os alunos)
void *desvioPadraoCidade(Regioes *r, double *dpCidades, int *medias) {
    int i,j;

    for(i = 0; i < reg->R * reg->C; i++){
        for(j = 0; i < reg->A; i++)
            DesvioPadrao[i] += (reg->m[i][j] - medias[i]) * (reg->m[i][j] - medias[i]);
        }
        DesvioPadrao[i] = sqrt(DesvioPadrao[i]/(double)(reg->R * reg->C));
    }
}

// Calcula a media aritmetica para cada regiao (entre as cidades)
void MediaAritmeticaRegiao(Regioes *reg, int *maRegiao) {
    int i, j, cont = 0;
    int regiao = 0;

    for (i = 0; i < reg->C*reg->R; i++) {
        for (j = 0; j < reg->A; j++) {
            maRegiao[regiao] += reg->m[i][j];
        }
        cont++;
        if ((cont%reg->C) == 0) {
            maRegiao[regiao] = maRegiao[regiao] / (double)reg->C*reg->A;
            regiao++;
        }
    }    
}

void exibe(Regioes *r) {
    printf("R = %d, C = %d, A = %d\n", r->R, r->C, r->A);
    int lin = r->R * r->C;
    int col = r->A;

    for (int i = 0; i < lin; i++) {
        for (int j = 0; j < col; j++) {
            printf("%d ", r->m[i][j]);
        }
        printf("\n");
    }
}

void libera_memoria(Regioes *r) {
    for (int i = 0; i < r->R * r->C; i++) {
        free(r->m[i]);
    }
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
    int *maioresCidade     = (int *)    calloc(r->R * r->C, sizeof(int));
    int *menoresCidade     = (int *)    calloc(r->R * r->C, sizeof(int));
    double *maCidade       = (double *) calloc(r->R * r->C, sizeof(double)); // maCidade = medias aritmeticas das cidades
    double *medianasCidade = (double *) calloc(r->R * r->C, sizeof(double));
    double *dpCidade       = (double *) calloc(r->R * r->C, sizeof(double)); // dpCidade = desvios padroes das cidades
    
    // Vetores que armazenam os dados para cada regiao
    double *maRegiao = (int *) malloc(sizeof(int)*reg->R);
    double maBrasil;
    
    
    // Ordena as notas de cada cidade (ordena as linhas)
    for(i = 0; i < regioes->C*regioes->R; i++) {
        quickSort(regioes->m[i], 0, regioes->A-1);
    }

    /* ..:: Chamando Funções ::.. */
    desvioPadraoCidade(Regioes *r, int *medias, double *dpCidades);
    MedianaCidade(Regioes *reg, double *medianasCidade);
    MediaAritmeticaCidade(Regioes *reg, double *maCidade);
    menorCidade(Regioes *r, int *menores);
    maiorCidade(Regioes *r, int *maiores);
    
    libera_memoria(regioes);
    return 0;
}