/*  Trabalho 1 - Programacao Concorrente - SSC0143
    Segundo semestre de 2019

    Daniel Penna Chaves Bertazzo - 10349561
    Gabriel Seiji Matsumoto      - 10295332
    Leonardo Miassi Netto        - 9326688
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

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

Regioes *le_teste() {
    int R, A, C;
    int i, j;
    scanf("%d %d %d", &R, &C, &A);

    Regioes *aux = (Regioes *) malloc(sizeof(Regioes));
    aux->R = R;
    aux->C = C;
    aux->A = A;

    aux->m = (int *) malloc(R*C*A * sizeof(int));
    for(i = 0; i < R*C; i++) {    
        for(j = 0; j < A; j++) {
            scanf("%d", &aux->m[i*A + j]);
        }
    }

    return aux;
}

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
            aux->m[i*A + j] = rand() % 101; // Valores de 0 a 100
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
        maioresCidade[i] = r->m[(i*r->A) + (r->A-1)];
    }
}


// Calcula as menores notas por cidade
void menorCidade(Regioes *r, int *menoresCidade) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    for (int i = 0; i < r->R*r->C; i++) {
        menoresCidade[i] = r->m[i*r->A];
    }
}


// Calcula media aritmetica para cada cidade (entre os alunos)
void mediaAritmeticaCidade(Regioes *r, double *maCidade) {
    int i,j;

    for (i = 0; i < r->C*r->R; i++) {
        for (j = 0; j < r->A; j++) {
            maCidade[i] += r->m[i*r->A + j];
        }
        maCidade[i] = maCidade[i] / (double)r->A;
    }
}

// Calcula mediana para cada cidade (entre os alunos)
void medianaCidade(Regioes *r, double *medianasCidade){
    int i,j;
    int decisao, mid;

    // meio = num_colunas / 2
    mid = r->A/2;

    if((r->A%2) == 0) decisao = 0;
    else decisao = 1;

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

    for(i = 0; i < r->R * r->C; i++) {
        for(j = 0; j < n; j++) {
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
    for (int i = 0; i < r->R; i++) {
        maioresRegiao[i] = r->m[(i*tamRegiao) + (tamRegiao-1)];
    }
}

// Calcula as menores notas por regiao
void menorRegiao(Regioes *r, int *menoresRegiao) {
    // Vetor (r->m[i]) ja esta ordenado => menor nota na primeira posicao
    int tamRegiao = r->C * r->A;
    for (int i = 0; i < r->R; i++) {
        menoresRegiao[i] = r->m[i*tamRegiao];
    }
}

// Calcula a media aritmetica para cada regiao (entre as cidades)
void mediaAritmeticaRegiao(Regioes *r, double *maRegiao) {
    int i, j;
    int tamRegiao = r->A * r->C;

    // Itera sobre as R regioes
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

    for(i = 0; i < r->R ; i++){
        for(j = 0; j < tamRegiao; j++) {
            dpRegiao[i] += (r->m[i*tamRegiao + j] - maRegiao[i]) * (r->m[i*tamRegiao + j] - maRegiao[i]);
        }
        dpRegiao[i] = sqrt(dpRegiao[i]/(double)(tamRegiao-1.0));
    }
}

/**********************************************************************************************************/
/* ..:: CALCULOS DO BRASIL ::.. */

int maiorBrasil(Regioes *r) {
    // Vetor (r->m) esta ordenado => maior na ultima posicao
    return r->m[(r->R * r->C * r->A) - 1];
}

int menorBrasil(Regioes *r) {
    // Vetor (r->m) esta ordenado => maior na primeira posicao
    return r->m[0];
}

double mediaAritmeticaBrasil(Regioes *r) {
    double maBrasil = 0.0;
    int tamBrasil = r->R * r->C * r->A; 

    for (int i = 0; i < tamBrasil; i++) {
        maBrasil += r->m[i];
    }

    return (maBrasil / (double)(tamBrasil));
}

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

double desvioPadraoBrasil(Regioes *r, double maBrasil) {
    double soma;
    double n = r->R*r->C*r->A;

    for (int i = 0; i < n; i++){
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

int main(int argc, char const *argv[]) {
    int i, j;
    clock_t inicio, fim;
    double delta;

    Regioes *regioes = le_entrada();
    // Regioes *regioes = le_teste();


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
    
    // Variaveis que armazenam os dados para o Brasil todo
    int majorBrasil = 0; 
    int minorBrasil = 0;
    double maBrasil = 0.0;
    double meBrasil = 0.0;
    double dpBrasil = 0.0;

    // printf("MATRIZ ORIGINAL: \n");
    // exibe(regioes);
    // printf("\n\n");

    inicio = clock();

    // Ordena as notas de cada cidade (ordena as linhas)
    int tamCidade = regioes->A;
    for (i = 0; i < regioes->C*regioes->R; i++) {
        quickSort(regioes->m, i * tamCidade, (i*tamCidade) + (tamCidade-1));
    }


    // printf("MATRIZ ORDENADA POR CIDADE: \n");
    // exibe(regioes);
    // printf("\n\n");

    // Chamando funcoes para as cidades
    maiorCidade(regioes, maioresCidade);
    menorCidade(regioes, menoresCidade);
    mediaAritmeticaCidade(regioes, maCidade);
    medianaCidade(regioes, medianasCidade);
    desvioPadraoCidade(regioes, dpCidade, maCidade);

    // Ordena as notas das regioes (ordena os blocos CxA que representam as regioes)
    int tamRegiao = regioes->C * regioes->A;
    for (i = 0; i < regioes->R; i++) {
        quickSort(regioes->m, i * tamRegiao, (i*tamRegiao) + (tamRegiao-1));
    }

    // printf("MATRIZ ORDENADA POR REGIAO: \n");
    // exibe(regioes);
    // printf("\n\n");

    // Chamando funcoes para as regioes
    maiorRegiao(regioes, maioresRegiao);
    menorRegiao(regioes, menoresRegiao);
    mediaAritmeticaRegiao(regioes, maRegiao);
    medianaRegiao(regioes, medianasRegiao);
    desvioPadraoRegiao(regioes, dpRegiao, maRegiao);

    // Ordena todos os dados do pais
    quickSort(regioes->m, 0, (regioes->R * regioes->C * regioes->A) - 1);

    // printf("MATRIZ ORDENADA POR Brasil: \n");
    // exibe(regioes);
    // printf("\n\n");


    // Chamando funcoes para o Brasil
    majorBrasil = maiorBrasil(regioes);
    minorBrasil = menorBrasil(regioes);
    maBrasil = mediaAritmeticaBrasil(regioes);
    meBrasil = medianaBrasil(regioes);
    dpBrasil = desvioPadraoBrasil(regioes,maBrasil);

    // Para de medir o tempo
    fim = clock();
    delta = ((double) (fim - inicio)) / CLOCKS_PER_SEC;

    /* ..:: Imprimir os resultados ::.. */
    double maiorMediaCidade = maCidade[0], maiorMediaRegiao = maRegiao[0];
    int qualCidade = 0, qualCidadeRegiao = 0, qualRegiao = 0;
    int n = regioes->A, cont = 0;
    
    for(i = 0; i < regioes->R; i++){
        for(j = 0; j < regioes->C; j++){
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

    //Resultado geral do Brasil
    printf("Brasil: menor: %d, maior: %d, mediana: %.2lf, media: %.2lf e DP: %.2lf\n",
            minorBrasil, majorBrasil, meBrasil, maBrasil, dpBrasil);

    printf("\n");
    //Melhores Classificações
    printf("Melhor região: Região %d\n", qualRegiao);
    printf("Melhor cidade: Região %d, Cidade %d\n", qualCidadeRegiao, qualCidade);


    printf("Tempo de resposta sem considerar E/S, em segundos: %fs\n", delta);

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