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
void mediaAritmeticaCidade(Regioes *r, double *maCidade) {
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
void medianaCidade(Regioes *r, double *medianasCidade){
    int i,j;
    int decisao, mid;

    mid = r->A/2;

    if((r->A%2) == 0) decisao = 0;
    else decisao = 1;

    for(i = 0; i < r->R * r->C; i++){
        if(decisao) {
            // medianasCidade[i] = r->m[i][mid+1];
            medianasCidade[i] = r->m[(i*r->A) + (mid+1)];
        }
        else {
            // medianasCidade[i] = (r->m[i][mid] + r->m[i][mid+1]) / 2.0;
            medianasCidade[i] = (r->m[i*r->A + mid] + r->m[(i+r->A) + (mid+1)]) / 2.0;
        }

    }
}

// Calcula o desvio padrao para cada cidade (entre os alunos)
void *desvioPadraoCidade(Regioes *r, double *dpCidade, double *maCidade) {
    int i, j;

    for(i = 0; i < r->R * r->C; i++) {
        for(j = 0; j < r->A; j++) {
            // dpCidade[i] += (r->m[i][j] - maCidade[i]) * (r->m[i][j] - maCidade[i]);
            dpCidade[i] += (r->m[i*r->A + j] - maCidade[i]) * (r->m[i*r->A + j] - maCidade[i]);
        }
        dpCidade[i] = sqrt(dpCidade[i]/(double)(r->R * r->C));
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
void mediaAritmeticaRegiao(Regioes *r, double *maRegiao) {
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
void medianaRegiao(Regioes *r, double *medianasRegiao){
    int i, j;
    int decisao, mid;

    mid = (r->A * r->C)/2;

    if (((r->A * r->C)%2) == 0) decisao = 0;
    else decisao = 1;

    for(i = 0; i < r->R; i++){
        if(decisao) {
            // medianasRegiao[i] = r->m[i][mid+1];
            medianasRegiao[i] = r->m[(i*r->A) + (mid+1)];
        }
        else { 
            // medianasRegiao[i] = (r->m[i][mid] + r->m[i][mid+1]) / 2.0;
            medianasRegiao[i] = (r->m[i*r->A + mid] + r->m[(i+r->A) + (mid+1)]) / 2.0;
        }

    }
}

// Calcula o desvio padrao para cada regiao 
void *desvioPadraoRegiao(Regioes *r, double *dpRegiao, double *maRegiao) {
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

double medianaBrasil(Regiao *r) {
    int tamBrasil = r->R * r->C * r->A;
    int mid = tamBrasil / 2;

    if (((r->R * r->C * r->A) % 2) == 0) {
        return (double)(r->m[mid] + r->m[mid+1]) / 2.0
    }
    else {
        return (double)r->m[mid];
    }

}

double dpBrasil(Regioes *r) {
    
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

void ImprimeResultado(Regioes *r, double *menoresCidade, double *maioresCidade, double *medianasCidade, double *maCidade, double *dpCidade, double *menoresRegiao, double *maioresRegiao, double *medianasRegiao, double *maRegiao, double *dpRegiao){
    float MaiorMediaCidade = maCidade[0], MaiorMediaRegiao = maRegiao[0];
    int QualCidade=0, QualCidadeRegiao =0, QualRegiao=0;
    
    //Resultado por Cidade
    for(int i = 0; i < r->R; i++){
        for(int j = 0; j < r->C; j++){
            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %d, media: %d e DP: %d\n", i, j, menoresCidade[j], maioresCidade[j], medianasCidade[j], maCidade[j], dpCidade[j]);
            if(MaiorMediaCidade < maCidade[j]){
                MaiorMediaCidade = maCidade[j];
                QualCidade = j;
                QualCidadeRegiao = i;
            }
        }
        printf("\n");
    }
    
    //Resultado por Regiao
    for(int i = 0; i < r->R; i++){
        printf("Reg %d: menor: %d, maior: %d, mediana: %d, media: %d e DP: %d\n", i, menoresRegiao[j], maioresRegiao[j], medianasRegiao[j], maRegiao[j], dpRegiao[j]);    
        if(MaiorMediaRegiao < maRegiao[i]){
            MaiorMediaRegiao = maRegiao[i];
            QualRegiao = i;
        }
    }

    //Resultado geral do Brasil
    
    //Melhores Classificações
    printf("Melhor região: Região %d\n", QualRegiao);
    printf("Melhor cidade: Região %d, Cidade %d\n", QualCidadeRegiao, QualCidade);

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
    
    // Variaveis que armazenam os dados para o Brasil todo
    int maiorBrasil      = 0; 
    int menorBrasil      = 0;
    double maBrasil      = 0.0 
    double medianaBrasil = 0.0;
    double dpBrasil      = 0.0;


    // Ordena as notas de cada cidade (ordena as linhas)
    int tamCidade = regioes->A;
    for (i = 0; i < regioes->C*regioes->R; i++) {
        quickSort(regioes->m, i * tamCidade, (i*tamCidade) + (tamCidade-1));
    }

    /* ..:: Chamando funcoes para as cidades ::.. */
    /* CHAMAR AS FUNCOES DAS CIDADES AQUI */    


    // Ordena as notas das regioes (ordena os blocos CxA que representam as regioes)
    int tamRegiao = regioes->C * regioes->A;
    for (i = 0; i < regioes->R; i++) {
        quickSort(regioes->m, i * tamRegiao, (i*tamRegiao) + (tamRegiao-1));
    }

    /* ..:: Chamando funcoes para as regioes ::.. */
    /* CHAMAR AS FUNCOES DAS REGIOES AQUI */    

    // Ordena todos os dados do pais
    quickSort(regioes->m, 0, regioes->R * regioes->C * regioes->A);

    /* ..:: Chamando funcoes para o Brasil ::.. */
    /* CHAMAR AS FUNCOES DO BRASIL AQUI */    


    /* ..:: Liberacao de memoria ::.. */
    libera_memoria(regioes);
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