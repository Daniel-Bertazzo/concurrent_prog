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
            printf("%d ", aux->m[i][j]);
        }
        printf("\n");
    }
    
    return aux;
}

void EncontraMedias(Regioes *reg){
    int i,j;
    int cont = 0;
    for(i = 0)
}

void libera_memoria(Regioes *r) {
    for (int i = 0; i < r->R * r->C; i++) {
        free(r->m[i]);
    }
    free(r->m);
    free(r);    
}

int main(int argc, char const *argv[]) {
    Regioes *regioes = le_entrada();

    libera_memoria(regioes);
    return 0;
}