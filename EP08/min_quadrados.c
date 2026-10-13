/*
 * Ajuste polinomial pelo Metodo dos Minimos Quadrados Discreto
 *
 * Entrada (teclado ou arquivo):
 *     n                -> numero de pontos (n <= 100)
 *     x1 f(x1)
 *     ...
 *     xn f(xn)
 *     m                -> grau do polinomio (m <= 3)
 *
 * Saida: coeficientes a0..am de P(x) = a0 + a1 x + ... + am x^m,
 *        erro (soma dos quadrados dos residuos) e tabela de residuos.
 *
 * Compilar: gcc minimos_quadrados.c -o mq -lm
 * Uso:      ./mq              (le do teclado)
 *           ./mq entrada.txt  (le do arquivo)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 100
#define MMAX 3

/* Resolve A*c = b (ordem k) por eliminacao de Gauss com pivotamento parcial.
   Retorna 0 em caso de sucesso, 1 se a matriz for singular. */
int resolve_sistema(int k, double A[MMAX + 1][MMAX + 2], double c[])
{
    int i, j, l, p;
    double aux, f;

    for (j = 0; j < k; j++) {
        p = j;
        for (i = j + 1; i < k; i++)
            if (fabs(A[i][j]) > fabs(A[p][j])) p = i;

        if (fabs(A[p][j]) < 1e-14) return 1;

        if (p != j)
            for (l = 0; l <= k; l++) {
                aux = A[j][l]; A[j][l] = A[p][l]; A[p][l] = aux;
            }

        for (i = j + 1; i < k; i++) {
            f = A[i][j] / A[j][j];
            for (l = j; l <= k; l++) A[i][l] -= f * A[j][l];
        }
    }

    for (i = k - 1; i >= 0; i--) {
        aux = A[i][k];
        for (j = i + 1; j < k; j++) aux -= A[i][j] * c[j];
        c[i] = aux / A[i][i];
    }
    return 0;
}

int main(int argc, char *argv[])
{
    FILE *in = stdin;
    int n, m, i, j, k, teclado = 1;
    double x[NMAX], y[NMAX];
    double S[2 * MMAX + 1];            /* S[p] = soma x^p      */
    double A[MMAX + 1][MMAX + 2];      /* matriz aumentada     */
    double a[MMAX + 1];
    double pot, soma_erro = 0.0, yc, r;

    if (argc > 1) {
        in = fopen(argv[1], "r");
        if (in == NULL) {
            printf("Erro: nao foi possivel abrir o arquivo '%s'.\n", argv[1]);
            return 1;
        }
        teclado = 0;
    }

    /* ---------- leitura de n ---------- */
    if (teclado) printf("Numero de pontos n (max %d): ", NMAX);
    if (fscanf(in, "%d", &n) != 1) { printf("Erro na leitura de n.\n"); return 1; }
    if (n < 1 || n > NMAX) {
        printf("Erro: n deve estar entre 1 e %d (n = %d).\n", NMAX, n);
        return 1;
    }

    /* ---------- leitura dos pontos ---------- */
    for (i = 0; i < n; i++) {
        if (teclado) printf("Ponto %d - x f(x): ", i + 1);
        if (fscanf(in, "%lf %lf", &x[i], &y[i]) != 2) {
            printf("Erro na leitura do ponto %d.\n", i + 1);
            return 1;
        }
    }

    /* ---------- leitura de m ---------- */
    if (teclado) printf("Grau do polinomio m (max %d): ", MMAX);
    if (fscanf(in, "%d", &m) != 1) { printf("Erro na leitura de m.\n"); return 1; }
    if (m < 0 || m > MMAX) {
        printf("Erro: m deve estar entre 0 e %d (m = %d).\n", MMAX, m);
        return 1;
    }
    if (n < m + 1) {
        printf("Erro: sao necessarios ao menos m+1 = %d pontos (n = %d).\n", m + 1, n);
        return 1;
    }
    if (!teclado) fclose(in);

    /* ---------- sistema normal ---------- */
    for (k = 0; k <= 2 * m; k++) {
        S[k] = 0.0;
        for (i = 0; i < n; i++) S[k] += pow(x[i], k);
    }
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= m; j++) A[i][j] = S[i + j];
        A[i][m + 1] = 0.0;
        for (j = 0; j < n; j++) A[i][m + 1] += pow(x[j], i) * y[j];
    }

    printf("\nSistema normal (matriz aumentada):\n");
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= m + 1; j++) printf("%16.8f ", A[i][j]);
        printf("\n");
    }

    if (resolve_sistema(m + 1, A, a)) {
        printf("\nErro: sistema singular (pontos insuficientes ou x repetidos "
               "para o grau pedido).\n");
        return 1;
    }

    /* ---------- resultados ---------- */
    printf("\nCoeficientes do polinomio de grau %d:\n", m);
    for (i = 0; i <= m; i++) printf("  a%d = %.10f\n", i, a[i]);

    printf("\nP(x) = ");
    for (i = 0; i <= m; i++) {
        printf("%+.6f", a[i]);
        if (i == 1) printf("*x");
        if (i > 1)  printf("*x^%d", i);
        printf(" ");
    }
    printf("\n");

    printf("\n%4s %14s %14s %14s %14s\n", "i", "x", "f(x)", "P(x)", "residuo");
    for (i = 0; i < n; i++) {
        yc = 0.0; pot = 1.0;
        for (j = 0; j <= m; j++) { yc += a[j] * pot; pot *= x[i]; }
        r = y[i] - yc;
        soma_erro += r * r;
        printf("%4d %14.6f %14.6f %14.6f %14.6e\n", i + 1, x[i], y[i], yc, r);
    }
    printf("\nErro (soma dos quadrados dos residuos) = %.10e\n", soma_erro);

    return 0;
}