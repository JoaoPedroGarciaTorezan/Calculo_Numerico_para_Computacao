/*Instituto de Matemática e Computação
  CMAC05 – Cálculo Numérico para Computação
  Exercício Prático 07 - 29/09/26
  João Pedro Garcia Torezan - 2025002063
  Rodrigo Silvestre Ribeiro de Oliveira - 2024001965

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100

static double px[MAX], py[MAX];
static int npontos, grau;

static void imprime_fator(double xi)
{
    if (xi >= 0) {
        printf("(x-%g)", xi);
    }
    else {
        printf("(x+%g)", -xi);
    }
}

static void preve(double x)
{
    double xs[MAX], tab[MAX][MAX], dist[MAX];
    double resumo_val[MAX], resumo_err[MAX];
    int i, j, k, g;

    /* Ordena os pontos pela distancia a x (mais proximo primeiro) */
    for (i = 0; i < npontos; i++) {
        double d = fabs(px[i] - x), yv = py[i], xv = px[i];
        j = i - 1;
        while (j >= 0 && dist[j] > d) {
            dist[j + 1] = dist[j];
            xs[j + 1]   = xs[j];
            tab[j + 1][0] = tab[j][0];
            j--;
        }
        dist[j + 1] = d;
        xs[j + 1]   = xv;
        tab[j + 1][0] = yv;
    }

    /* Diferencas divididas na ordem da sequencia escolhida */
    for (j = 1; j < npontos; j++) {
        for (i = 0; i < npontos - j; i++) {
            tab[i][j] = (tab[i + 1][j - 1] - tab[i][j - 1]) / (xs[i + j] - xs[i]);
        }
    }

    printf("\n================ Previsao para x = %g ================\n", x);
    printf("Sequencia de pontos (do mais proximo ao mais distante):\n");
    for (i = 0; i < npontos; i++) {
        printf("  x%d = %g   f = %g\n", i, xs[i], tab[i][0]);
    }

    printf("\nTabela de diferencas divididas:\n");
    for (i = 0; i < npontos; i++) {
        printf("x%-2d = %8.4f | ", i, xs[i]);
        for (j = 0; j < npontos - i; j++) {
            printf("%12.6f ", tab[i][j]);
        }
        printf("\n");
    }

    for (g = 1; g <= grau; g++) {
        double val, produto = 1.0;

        printf("\n--- Polinomio de grau %d (pontos x0..x%d) ---\n", g, g);
        printf("n(x) = %g", tab[0][0]);
        for (k = 1; k <= g; k++) {
            printf(" + %g.", tab[0][k]);
            for (i = 0; i < k; i++) {
                imprime_fator(xs[i]);
            }
        }
        printf("\n");

        /* Horner */
        val = tab[0][g];
        for (k = g - 1; k >= 0; k--) {
            val = val * (x - xs[k]) + tab[0][k];
        }
        printf("n%d(%g) = %.6f\n", g, x, val);
        resumo_val[g] = val;

        /* Estimativa de erro desta previsao */
        double dnext = tab[0][g + 1], erro;
        for (i = 0; i <= g; i++) {
            produto *= (x - xs[i]);
        }
        erro = fabs(produto * dnext);
        printf("Proximo ponto da tabela: x%d = %g\n", g + 1, xs[g + 1]);
        printf("|(x-x0)...(x-x%d)| = %.6e ; |D%d| = %.6e\n", g, fabs(produto), g + 1, fabs(dnext));
        printf("|E%d(%g)| ~= %.6e\n", g, x, erro);
        resumo_err[g] = erro;
    }

    printf("\n----------- Resumo para x = %g -----------\n", x);
    printf("Grau |      n(x)      | Erro estimado\n");
    for (g = 1; g <= grau; g++) {
        printf("%4d | %14.6f | %.6e\n", g, resumo_val[g], resumo_err[g]);
    }
}

int main(void)
{
    int i;
    double x;
    char op;

    printf("Quantidade de pontos: ");
    if (scanf("%d", &npontos) != 1 || npontos < 2 || npontos > MAX) {
        printf("Quantidade de pontos invalida (2 a %d).\n", MAX);
        return 1;
    }

    printf("Grau maximo desejado (polinomios de grau 1 ate este): ");
    if (scanf("%d", &grau) != 1 || grau < 1 || grau > npontos - 2) {
        printf("Grau invalido: deve estar entre 1 e %d (o erro usa o proximo ponto).\n", npontos - 2);
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
    for (i = 0; i < npontos; i++) {
        int j;
        for (j = i + 1; j < npontos; j++)
            if (px[i] == px[j]) {
                printf("Erro: existem valores de x repetidos.\n");
                return 1;
            }
    }

    do {
        printf("\nValor de x a interpolar: ");
        if (scanf("%lf", &x) != 1) {
            printf("Valor de x invalido.\n");
            return 1;
        }
        preve(x);
        printf("\nDeseja fazer outra previsao? (s/n): ");
        if (scanf(" %c", &op) != 1) {
            break;
        }
    } while (op == 's' || op == 'S');

    return 0;
}