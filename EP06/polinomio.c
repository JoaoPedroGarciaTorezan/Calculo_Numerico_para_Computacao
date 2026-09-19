/*Instituto de Matemática e Computação
CMAC05 – Cálculo Numérico para Computação
Exercício Prático 06 - 15/09/26 
João Pedro Garcia Torezan - 2025002063
João Otávio Zampieri Fornazeiro - 2025005163
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PONTOS 10

/* Resolve o sistema linear A * a = b usando Eliminação de Gauss com pivotamento parcial. */
int resolverSistema(int n, double A[MAX_PONTOS+1][MAX_PONTOS+2], double a[MAX_PONTOS+1]) {
    int i, j, k, piv;
    int tam = n + 1;

    for (k = 0; k < tam; k++) {
        piv = k;
        for (i = k + 1; i < tam; i++) {
            if (fabs(A[i][k]) > fabs(A[piv][k]))
                piv = i;
        }
        if (fabs(A[piv][k]) < 1e-12) {
            printf("Sistema singular ou pontos com x repetidos!\n");
            return 0;
        }
        // Troca linhas se necessário 
        if (piv != k) {
            for (j = 0; j <= tam; j++) {
                double tmp = A[k][j];
                A[k][j] = A[piv][j];
                A[piv][j] = tmp;
            }
        }

        // Eliminação 
        for (i = k + 1; i < tam; i++) {
            double fator = A[i][k] / A[k][k];
            for (j = k; j <= tam; j++)
                A[i][j] -= fator * A[k][j];
        }
    }

    //Substituição regressiva 
    for (i = tam - 1; i >= 0; i--) {
        double soma = A[i][tam]; //  termo b 
        for (j = i + 1; j < tam; j++)
            soma -= A[i][j] * a[j];
        a[i] = soma / A[i][i];
    }
    return 1;
}

// Avalia o polinômio P(x) = a0 + a1*x + ... + an*x^n usando Horner
double avaliarPolinomio(double a[], int n, double x) {
    double resultado = a[n];
    int i;
    for (i = n - 1; i >= 0; i--) {
        resultado = resultado * x + a[i];
    }
    return resultado;
}

int main() {
    double x[MAX_PONTOS], y[MAX_PONTOS];
    double A[MAX_PONTOS+1][MAX_PONTOS+2]; // Matriz aumentada
    double a[MAX_PONTOS+1]; //Matriz de coeficientes
    int n, i, j, grau;
    double xAvaliar, xmin, xmax, resultado;

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

    grau = n - 1;

    // Monta a matriz de Vandermonde
    for (i = 0; i < n; i++) {
        double potencia = 1.0;
        for (j = 0; j <= grau; j++) {
            A[i][j] = potencia;
            potencia *= x[i];
        }
        A[i][grau + 1] = y[i]; // termo independente b
    }

    if (!resolverSistema(grau, A, a)) {
        printf("Nao foi possivel calcular os coeficientes.\n");
        return 1;
    }

    printf("\nPolinomio interpolador de grau %d:\n", grau);
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

    // Determina o intervalo
    xmin = x[0];
    xmax = x[0];
    for (i = 1; i < n; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
    }

    // Solicita o ponto para avaliação, validando o intervalo 
    do {
        printf("\nDigite o valor de x para calcular P%d(x) (intervalo: %.4f a %.4f): ",
               grau, xmin, xmax);
        scanf("%lf", &xAvaliar);
        if (xAvaliar < xmin || xAvaliar > xmax)
            printf("Valor fora do intervalo! Tente novamente.\n");
    } while (xAvaliar < xmin || xAvaliar > xmax);

    resultado = avaliarPolinomio(a, grau, xAvaliar);
    printf("\nP%d(%.4f) = %.6f\n", grau, xAvaliar, resultado);

    return 0;
}