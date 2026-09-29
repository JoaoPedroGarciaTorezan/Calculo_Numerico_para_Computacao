#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100

/* Imprime o fator (x - xi) tratando o sinal de xi */
static void imprime_fator(double xi)
{
    if (xi >= 0)
        printf("(x-%g)", xi);
    else
        printf("(x+%g)", -xi);
}

int main(void)
{
    int npontos, grau, m, i, j, k, inicio;
    double x, px[MAX], py[MAX];
    double xs[MAX], tab[MAX][MAX];

    printf("Quantidade de pontos: ");
    if (scanf("%d", &npontos) != 1 || npontos < 1 || npontos > MAX) {
        printf("Quantidade de pontos invalida (1 a %d).\n", MAX);
        return 1;
    }

    printf("Grau do polinomio desejado: ");
    if (scanf("%d", &grau) != 1 || grau < 0 || grau > npontos - 1) {
        printf("Grau invalido: deve estar entre 0 e %d.\n", npontos - 1);
        return 1;
    }

    printf("Valor de x a interpolar: ");
    if (scanf("%lf", &x) != 1) {
        printf("Valor de x invalido.\n");
        return 1;
    }

    printf("\nInforme os pontos (xi yi):\n");
    for (i = 0; i < npontos; i++) {
        printf("Ponto %d: ", i);
        if (scanf("%lf %lf", &px[i], &py[i]) != 2) {
            printf("Entrada invalida.\n");
            return 1;
        }
    }

    // Ordena os pontos por x 
    for (i = 1; i < npontos; i++) {
        double tx = px[i], ty = py[i];
        j = i - 1;
        while (j >= 0 && px[j] > tx) {
            px[j + 1] = px[j];
            py[j + 1] = py[j];
            j--;
        }
        px[j + 1] = tx;
        py[j + 1] = ty;
    }

    // Verifica x repetidos (divisao por zero nas diferencas divididas) 
    for (i = 1; i < npontos; i++) {
        if (px[i] == px[i - 1]) {
            printf("Erro: existem valores de x repetidos.\n");
            return 1;
        }
    }

    // Seleciona grau+1 pontos consecutivos mais proximos de x 
    m = grau + 1;
    inicio = 0;
    {
        double melhor = INFINITY;
        for (i = 0; i + m <= npontos; i++) {
            double c1 = fabs(x - px[i]);
            double c2 = fabs(x - px[i + m - 1]);
            double custo = c1 > c2 ? c1 : c2;
            if (custo < melhor) {
                melhor = custo;
                inicio = i;
            }
        }
    }

    for (i = 0; i < m; i++) {
        xs[i] = px[inicio + i];
        tab[i][0] = py[inicio + i];
    }

    // Tabela de diferencas divididas 
    for (j = 1; j < m; j++)
        for (i = 0; i < m - j; i++)
            tab[i][j] = (tab[i + 1][j - 1] - tab[i][j - 1]) / (xs[i + j] - xs[i]);

    // Impressao da tabela
    printf("\nPontos utilizados e tabela de diferencas divididas:\n");
    for (i = 0; i < m; i++) {
        printf("x%d = %10.4f | ", i, xs[i]);
        for (j = 0; j < m - i; j++)
            printf("%12.6f ", tab[i][j]);
        printf("\n");
    }

    // Coeficientes Dk = tab[0][k] 
    printf("\nCoeficientes:\n");
    for (k = 0; k < m; k++)
        printf("D%d = %.6f\n", k, tab[0][k]);

    // Impressao do polinomio 
    printf("\nn(x) = %g", tab[0][0]);
    for (k = 1; k < m; k++) {
        printf(" + %g.", tab[0][k]);
        for (i = 0; i < k; i++)
            imprime_fator(xs[i]);
    }
    printf("\n");

    // Avaliacao em x (forma de Horner) 
    {
        double resultado = tab[0][m - 1];
        for (k = m - 2; k >= 0; k--)
            resultado = resultado * (x - xs[k]) + tab[0][k];
        printf("\nn(%g) = %.6f\n", x, resultado);
    }

    return 0;
}