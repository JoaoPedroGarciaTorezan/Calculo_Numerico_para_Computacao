#include <stdio.h>
#include <math.h>

double f(double x) {
    double ln = log(x);
    return x*ln - 2;
}

int main () {

    double a, b, e, x;

    printf("Limite inferior do intervalo da raiz: ");
    scanf("%lf", &a);
    printf("Limite superior do intervalo da raiz: ");
    scanf("%lf", &b);
    printf("Prescisao desejada: ");
    scanf("%lf", &e);

    int k = 1;
    while((b - a) >= e){

        double m= f(a);
        x = (a+b)/2;

        if(m*f(x) > 0)
            a = x;
        else b = x;
        k++;
    }

    x = (x+b)/2;

    printf("Raiz aproximada = %.4lf\n", x);
    printf("Número de iterações: %d\n", k);

    return 0;
}