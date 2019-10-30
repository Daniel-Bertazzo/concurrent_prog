/*  Trabalho 1 - Programacao Concorrente - SSC0143
    Segundo semestre de 2019

    Daniel Penna Chaves Bertazzo - 10349561
    Gabriel Seiji Matsumoto      - 10295332
    Leonardo Miassi Netto        - 

*/

#include <stdlib.h>
#include <stdio.h>

// Struct que representa as regioes
typedef struct Regioes {
    int R, C, A; // Numero de regioes, cidades e alunos por cidade
    int **m; // Matriz que representa os dados da regiao
} Regioes;


/* ..:: QUICKSORT RETIRADO DO GEEKS FOR GEEKS ::.. */
/********************************************************************************************/

void swap(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
} 
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (int arr[], int low, int high) 
{ 
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
void quickSort(int arr[], int low, int high) 
{ 
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

int *maior(Regioes *r) {
    // Armazena as maiores notas de cada cidade
    int *vet = (int *) malloc(sizeof(int) * r->R * r->C);

    // Vetor ja esta ordenado => maior nota na ultima posicao
    for (int i = 0; i < r->R*r->C; i++) {
        vet[i] = r->m[i][A-1];
    }
    
    return vet;
}

int *menor(Regioes *r){
    // Vetor para armazenar as menores notas de cada cidade
    int *vet = (int *)malloc(sizeof(int)*r->R*r->C);
    
    // Vetor ja esta ordenado => menor nota na primeira posicao
    for (int i = 0; i < r->R*r->C; i++) {
        vet[i] = r->m[i][0];
    }

    return vet;
    
}

// Calcula media aritmetica para cada cidade (entre os alunos)
int *MediaAritmetica(Regioes *reg) {
    int i,j;
    // vet -> armazena as medias de cada linha da matriz (cada escola)
    int *vet = (int *)malloc(sizeof(int)*reg->C*reg->R);
    /* soma -> armazena a soma dos valores das linhas (notas dos alunos),
       para realizar o calculo da media */
    int soma = 0;

    for (i = 0; i < reg->C*reg->R; i++) {
        for (j = 0; j < reg->A; j++) {
            soma += reg->m[i][j];
        }
        vet[i] = soma / reg->A;
        soma = 0;
    }

    return vet;
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
    int i,j;
    Regioes *regioes = le_entrada();
    for(int j = 0; j < regioes->C*regioes->R; j++)
        quickSort(regioes->m[j], 0, regioes->A-1);


    libera_memoria(regioes);
    return 0;
}