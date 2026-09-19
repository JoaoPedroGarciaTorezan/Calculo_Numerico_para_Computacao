/*/*Instituto de Matemática e Computação
CMAC05 – Cálculo Numérico para Computação
Exercício Prático 06 - 15/09/26 
João Pedro Garcia Torezan - 2025002063
João Otávio Zampieri Fornazeiro - 2025005163
*/

#include <stdio.h>
#include <math.h>

#define MAX_PONTOS 10

/* Multiplica o polinômio 'poli' (grau atual = grauAtual, coeficientes
   em poli[0..grauAtual]) pelo binômio (x - r).
   Resultado sobrescreve 'poli' e o grau aumenta em 1. */
void multiplicarPorBinomio(double poli[], int *grauAtual, double r) {
    int i;
    /* nova posição de maior grau recebe o coeficiente antigo de maior grau */
    poli[*grauAtual + 1] = poli[*grauAtual];

    /* demais posições: poli_novo[i] = poli_antigo[i-1] - r * poli_antigo[i] */
    for (i = *grauAtual; i >= 1; i--) {
        poli[i] = poli[i - 1] - r * poli[i];
    }
    poli[0] = -r * poli[0];

    (*grauAtual)++;
}

/* Calcula P(x) usando Horner a partir dos coeficientes já expandidos */
double avaliarPolinomio(double a[], int n, double x) {
    double resultado = a[n];
    int i;
    for (i = n - 1; i >= 0; i--)
        resultado = resultado * x + a[i];
    return resultado;
}

int main() {
    double x[MAX_PONTOS], y[MAX_PONTOS];
    double a[MAX_PONTOS];          /* coeficientes finais do polinômio interpolador */
    double Li[MAX_PONTOS];         /* coeficientes do Li(x) expandido (numerador) */
    int n, i, j, k;
    double xAvaliar, xmin, xmax;

    printf("=== Interpolacao Polinomial - Formula de Lagrange ===\n\n");

    do {
        printf("Digite o numero de pontos (2 a %d): ", MAX_PONTOS);
        scanf("%d", &n);
        if (n < 2 || n > MAX_PONTOS)
            printf("Valor invalido! Informe entre 2 e %d.\n", MAX_PONTOS);
    } while (n < 2 || n > MAX_PONTOS);

    printf("\nDigite os %d pontos (x y):\n", n);
    for (i = 0; i < n; i++) {
        printf("Ponto %d: ", i + 1);
        scanf("%lf %lf", &x[i], &y[i]);
    }

    /* Verifica se ha valores de x repetidos */
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(x[i] - x[j]) < 1e-12) {
                printf("Erro: valores de x repetidos (x%d = x%d). Encerrando.\n", i, j);
                return 1;
            }
        }
    }

    int grau = n - 1;

    /* Zera o vetor de coeficientes finais */
    for (i = 0; i <= grau; i++)
        a[i] = 0.0;

    /* --- Constrói cada Li(x) expandido e acumula em a[] --- */
    for (i = 0; i < n; i++) {
        int grauLi = 0;
        double denominador = 1.0;

        Li[0] = 1.0; /* começa como polinômio constante = 1 */

        for (k = 0; k <= grau; k++) Li[k] = 0.0;
        Li[0] = 1.0;

        for (j = 0; j < n; j++) {
            if (j != i) {
                /* multiplica Li(x) por (x - x[j]) */
                multiplicarPorBinomio(Li, &grauLi, x[j]);
                /* acumula o denominador: (x[i] - x[j]) */
                denominador *= (x[i] - x[j]);
            }
        }

        /* soma y[i]*Li(x) / denominador aos coeficientes finais */
        for (k = 0; k <= grau; k++) {
            a[k] += (y[i] / denominador) * Li[k];
        }
    }

    /* Exibe o polinômio */
    printf("\nPolinomio interpolador de grau %d (via Lagrange):\n", grau);
    printf("P%d(x) = ", grau);
    for (i = grau; i >= 0; i--) {
        if (i == 0)
            printf("%+.6f", a[i]);
        else if (i == 1)
            printf("%+.6f*x ", a[i]);
        else
            printf("%+.6f*x^%d ", a[i], i);
    }
    printf("\n\nCoeficientes (a0 ... a%d):\n", grau);
    for (i = 0; i <= grau; i++)
        printf("a%d = %.6f\n", i, a[i]);

    /* Determina xmin e xmax entre os pontos fornecidos */
    xmin = x[0];
    xmax = x[0];
    for (i = 1; i < n; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
    }

    /* Solicita o ponto para avaliação, validando o intervalo */
    do {
        printf("\nDigite o valor de x para calcular P%d(x) (intervalo: %.4f a %.4f): ",
               grau, xmin, xmax);
        scanf("%lf", &xAvaliar);
        if (xAvaliar < xmin || xAvaliar > xmax)
            printf("Valor fora do intervalo! Tente novamente.\n");
    } while (xAvaliar < xmin || xAvaliar > xmax);

    double resultado = avaliarPolinomio(a, grau, xAvaliar);
    printf("\nP%d(%.4f) = %.6f\n", grau, xAvaliar, resultado);

    return 0;
}