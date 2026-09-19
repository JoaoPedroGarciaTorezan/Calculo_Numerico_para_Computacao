#include <stdio.h>
#include <math.h>

double f(double x) {
    double ln = log(x);
    return x*ln - 2;
}

int main () {

    double a, b, e1, e2, x;
    int k;

    printf("Limite inferior do intervalo da raiz: ");
    scanf("%lf", &a);
    printf("Limite superior do intervalo da raiz: ");
    scanf("%lf", &b);
    printf("Prescisao desejada para (b-a): ");
    scanf("%lf", &e1);
    printf("Prescisão desejada para f(x): ");
    scanf("%lf", &e2);

    if(b-a < e1){

        if(fabs(f(a)) < e2 || fabs(f(b)) < e2){

            if(fabs(f(a)) < e2)
                x = a;
            else x = 1;
        }
    }
    else{
        k=1;
        x = (a * f(b) - b * f(a)) / (f(b) - f(a));
        while((b-a) >= e1 || fabs(f(x)) >= e2) {

            if(f(a)*f(b) > 0)
                a = x;
            else b = x;

            x = (a * f(b) - b * f(a)) / (f(b) - f(a));
            k++;
        }
    }

    printf("Raiz aproximada = %.4lf\n", x);
    printf("Número de iterações: %d\n", k);
    return 0;
}
