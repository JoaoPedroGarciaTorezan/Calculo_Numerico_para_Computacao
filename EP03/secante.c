#include <stdio.h>
#include <math.h>
/*Instituto de Matemática e Computação
CMAC05 – Cálculo Numérico para Computação
Exercício Prático 03 - 25/08/26
Caio Sano Veiga - 2025002733
João Pedro Garcia Torezan - 2025002063
 */

/* Exercício 1
a) x * x - tan(x)
b) 2*cos(x) - exp(2*x)
c) x*x*x*x*x - 6 
Exercício 2
x*x*x - 2 *x*x - 3*x + 10*/
double f(double x) {
    // Função a ser escolhida
    return x * x - tan(x);
}

// Função que calcula a derivada numérica em um ponto x
double df(double x) {
    double h = 1e-5;
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

int main() {

    double x2, x1, x0, e1, e2;
    int k;

    printf("=========== Algoritmo da Secante ===========\n");
    printf("Valor do inicio do intervalo da raiz(x0): ");
    scanf("%lf", &x0);

    printf("Valor do final do intervalo da raiz(x1): ");
    scanf("%lf", &x1);

    printf("Precisao desejada para (b-a) (e1): ");
    scanf("%lf", &e1);

    printf("Precisao desejada para f(x) (e2): ");
    scanf("%lf", &e2);


    if(fabs(f(x0)) < e2)
        x2 = x0;
    else {
        if(fabs(f(x1)) < e2 || fabs(x1 - x0) < e1) {
            x2 = x1;
        }
        else {
            k = 1;
            x2 = x1 - (f(x1)/(f(x1) - f(x0))) * (x1 - x0);
            do {
                x0 = x1;
                k++;
                x1 = x2;
                x2 = x1 - (f(x1)/(f(x1) - f(x0))) * (x1 - x0);
            } while (fabs(x2 - x1) >= e1 && fabs(f(x2)) >= e2);
        }
    }
    printf("Proxima estimativa da raiz(x2): %lf\n", x2);
    printf("Quantidade de iteracoes: %d\n", k);

    return 0;
}