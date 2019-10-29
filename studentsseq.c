#include <stdlib.h>
#include <stdio.h>

// Struct que representa as regioes
typedef struct Regiao {
    int R, C, A; // Numero de cidades e alunos por cidade
    int **m; // Matriz que representa os dados da regiao
} Regiao;

Regiao *le_entrada() {
    // Le a entrada do arquivo
    int R, A, C, seed;
    int i, j;
    scanf("%d %d %d %d", &R, &C, &A, &seed);

    // Aloca as regioes e atribui os valores
    Regiao *aux = (Regiao *) malloc(sizeof(Regiao));
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
}

int main(int argc, char const *argv[]) {
    Regiao *regioes = le_entrada();

    //void libera_memoria()
    return 0;
}