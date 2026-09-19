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
x*x*x - 2 *x*x - 3*x + 10
Exercício 3
-2*sin(x) + x*/
double f(double x) {
    // Função a ser escolhida
    return x*x*x - 3*x + 1;
}

// Função que calcula a derivada numérica em um ponto x
double df(double x) {
    double h = 1e-5;
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

int main() {

    double x1, x0, e1, e2;
    int k;

    printf("=========== Método de Newton-Raphson ===========\n");
    printf("Estimativa inicial da raiz(x0): ");
    scanf("%lf", &x0);

    printf("Precisao desejada para (b-a) (e1): ");
    scanf("%lf", &e1);

    printf("Precisao desejada para f(x) (e2): ");
    scanf("%lf", &e2);

    if(fabs(df(x0)) > 1e-12)  {
        //É continua
        if(fabs(f(x0)) < e2){
            //se sim
            x1 = x0;
        }
        else{
            k=1;
            x1 = x0 - f(x0) / df(x0);
            do{
                x0 = x1;
                k++;
                x1 = x0 - f(x0) / df(x0);
            } while(fabs(x1 - x0) >= e1 && fabs(f(x1)) >= e2);
        }
        printf("Proxima estimativa da raiz(x1): %lf\n", x1);
        printf("Quantidade de iteracoes: %d\n", k);
    }
    else {
        printf("Nao pode ser calculado pelo metodo. \n");
    }


    return 0;
}