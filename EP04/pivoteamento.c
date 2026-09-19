// Instituto de Matemática e Computação
// CMAC05 – Cálculo Numérico para Computação
// Exercício Prático 04 - 01/09/26
// João Pedro Garcia Torezan - 2025002063
// Davi da Silva Lorena - 2025004264
#include <stdio.h>
#include <math.h>

#define N 10

void imprime(double a[N][N+1], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%8.3f ", a[i][j]);
        }
        printf("\n");
    }
}

void troca(double a[N][N+1], int n, int i, int j) {
    for (int k = 0; k <= n; k++) {
        double t = a[i][k]; a[i][k] = a[j][k]; a[j][k] = t;
    }
}

int main() {
    int n;
    double a[N][N+1], x[N], orig[N][N+1];
    int correta = 1;
    double tol = 1e-6;
    double soma, m, residuo;

    printf("Digite o valor de n: "); 
    scanf("%d", &n);
    printf("Digite a matriz aumentada [A|b]:\n");
    for (int i = 0; i < n; i++) {
        printf("Valores da equacao %d: ", i+1);
        for (int j = 0; j <= n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    // guarda a matriz original (antes do pivoteamento/eliminacao) para a verificacao final
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            orig[i][j] = a[i][j];

    printf("\nMatriz inicial:\n");
    imprime(a, n);

    for (int k = 0; k < n - 1; k++) {
        // busca maior pivo na coluna k
        int piv = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i][k]) > fabs(a[piv][k])) 
                piv = i;
        }

        if (fabs(a[piv][k]) < 1e-12) {
            printf("\nColuna %d nula. Sistema singular.\n", k+1);
            return 1;
        }

        if (piv != k) {
            printf("\nTroca linha %d <-> linha %d\n", k+1, piv+1);
            troca(a, n, k, piv);
        }

        for (int i = k + 1; i < n; i++) {
            m = a[i][k] / a[k][k];
            printf("linha%d = linha%d - (%.3f)*linha%d\n", i+1, i+1, m, k+1);
            for (int j = k; j <= n; j++) {
                a[i][j] -= m * a[k][j];
            }
        }
        printf("Matriz apos etapa %d:\n", k+1);
        imprime(a, n);
    }

    if (fabs(a[n-1][n-1]) < 1e-12) {
        printf("\nSistema sem solucao unica.\n");
        return 1;
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }

    printf("\nSolucao:\n");
    for (int i = 0; i < n; i++) {
        if(i==0)
            printf("x = %.4f\n", x[i]);
        if(i==1)
            printf("y = %.4f\n", x[i]);
        if(i==2)
            printf("z = %.4f\n", x[i]);
        if(i==3)
            printf("w = %.4f\n", x[i]);
    }

    //Verificacao
    printf("\nVerificaca:\n");
    for (int i = 0; i < n; i++) {
        soma = 0.0;
        for (int j = 0; j < n; j++) {
            soma += orig[i][j] * x[j];
        }
        residuo = soma - orig[i][n];
        printf("Equacao %d: soma calculada = %.6f | valor esperado = %.6f \n",
               i+1, soma, orig[i][n]);
        if (fabs(residuo) > tol) {
            correta = 0;
        }
    }

    if (correta) {
        printf("\nResultado verificado: a solucao satisfaz o sistema original.\n");
    } else {
        printf("\nAtencao: a solucao NAO satisfaz o sistema original dentro da tolerancia (%.0e).\n", tol);
    }

    return 0;
}