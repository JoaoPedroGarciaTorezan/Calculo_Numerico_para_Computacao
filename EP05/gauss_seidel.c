/*
Diego Pereira Nadur Castro - 2025003928
João Pedro Garcia Torezan - 2025002063
*/
#include <stdio.h>
#include <math.h>

int main() {
    int n, max_iter, i, j, k, p;
    double eps, max_delta, max_x, delta_r;
    double a[10][11];
    double x[10], x_antes[10];
    int houve_troca = 0;

    printf("Digite a quantidade de equacoes: ");
    scanf("%d", &n);

    printf("\nDigite os elementos da matriz:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    printf("\nDigite a aproximacao inicial (x0):\n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &x[i]);
    }

    printf("\nDigite o criterio de parada (eps): ");
    scanf("%lf", &eps);

    printf("Digite o numero maximo de iteracoes: ");
    scanf("%d", &max_iter);

    printf("\n=== SISTEMA INICIAL INSERIDO ===\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.4lf ", a[i][j]);
        }
        printf(" | %8.4lf\n", a[i][n]);
    }

     
    // Reordenamento automatico das linhas
    for (i = 0; i < n; i++) {
        int linha_maior = i;
        double maior_valor = fabs(a[i][i]);

        for (p = i + 1; p < n; p++) {
            if (fabs(a[p][i]) > maior_valor) {
                maior_valor = fabs(a[p][i]);
                linha_maior = p;
            }
        }

        if (linha_maior != i) {
            houve_troca = 1;
            printf("\n Troca automatica: Linha %d <-> Linha %d ", i + 1, linha_maior + 1);
            for (j = 0; j <= n; j++) {
                double temp = a[i][j];
                a[i][j] = a[linha_maior][j];
                a[linha_maior][j] = temp;
            }
        }
    }

    //Se o sistema foi modificado, imprime a nova matriz 
    if (houve_troca) {
        printf("\n Sistema Equivalente: \n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%8.4lf ", a[i][j]);
            }
            printf(" | %8.4lf\n", a[i][n]);
        }
    }

    // Iteracoes de Gauss-Seidel 
    for (k = 1; k <= max_iter; k++) {
        printf("\nIteracao %d \n", k);
        max_delta = 0.0;
        max_x = 0.0;

        for (i = 0; i < n; i++) {
            x_antes[i] = x[i];
        }

        for (i = 0; i < n; i++) {
            double soma = 0.0;
            printf("x%d = (%.4lf", i + 1, a[i][n]);

            for (j = 0; j < n; j++) {
                if (j != i) {
                    soma = soma + a[i][j] * x[j];
                    printf(" - (%.4lf * %.4lf)", a[i][j], x[j]);
                }
            }
            x[i] = (a[i][n] - soma) / a[i][i];
            printf(") / %.4lf = %.4lf\n", a[i][i], x[i]);
        }

        printf("\nErros absolutos:\n");
        for (i = 0; i < n; i++) {
            double delta = fabs(x[i] - x_antes[i]);
            printf("Delta x%d = %.4lf\n", i + 1, delta);
            if (delta > max_delta) {
                max_delta = delta;
            }
            if (fabs(x[i]) > max_x) {
                max_x = fabs(x[i]);
            }
        }

        delta_r = max_delta / max_x;
        printf("Erro relativo (Delta r) = %.4lf / %.4lf = %.4lf\n", max_delta, max_x, delta_r);

        if (delta_r < eps) {
            printf("\nCriterio de parada atingido na iteracao %d (Delta r < %.4lf).\n", k, eps);
            break;
        }
    }


    printf("\nResposta: \n");
    for (i = 0; i < n; i++) {
        printf("x%d = %.4lf\n", i + 1, x[i]);
    }

    return 0;
}