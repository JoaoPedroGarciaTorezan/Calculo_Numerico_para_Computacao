/*Instituto de Matemática e Computação
CMAC05 – Cálculo Numérico para Computação
Exercício Prático 10 - 06/11/2026
João Pedro Garcia Torezan - 2025002063
*/

#include <stdio.h>
#include <math.h>
#define PONTOS_VARREDURA 100000

// 1) Função a ser integrada f(x) 
double f(double x)
{
    return x*x*x*x;                
}

// 2) Segunda derivada f''(x) (necessária para o cálculo do erro)
double f2(double x)
{
    return 12*x*x;                 
}

/* 3) (OPCIONAL) Primitiva F(x), para comparar com o valor exato.
 *    Se não souber a primitiva, defina CONHECE_PRIMITIVA como 0.     */
#define CONHECE_PRIMITIVA 1        /* <-- ALTERAR: 1 = sim, 0 = não */

double F(double x)
{
    return exp(x);             
}

/* Trapézios simples (n = 1) */
double trapezio_simples(double a, double b)
{
    double h = b - a;
    return (h / 2.0) * (f(a) + f(b));
}

/* Trapézios generalizado (n > 1) */
double trapezio_generalizado(double a, double b, int n)
{
    double h = (b - a) / n;
    double soma = f(a) + f(b);
    int i;
    for (i = 1; i < n; i++) {
        soma += 2.0 * f(a + i * h);
    }
    return (h / 2.0) * soma;
}

/* Estima max |f''(x)| em [a,b] por varredura */
double max_f2(double a, double b)
{
    double max = fabs(f2(a)), v, x;
    int i;
    for (i = 1; i <= PONTOS_VARREDURA; i++) {
        x = a + (b - a) * i / PONTOS_VARREDURA;
        v = fabs(f2(x));
        if (v > max) {
            max = v;
        }
    }
    return max;
}

int main(void)
{
    double a, b, h, I, M, erro;
    int n;

    printf("=== Integracao pelo Metodo dos Trapezios ===\n");
    printf("Limite inferior a: ");
    if (scanf("%lf", &a) != 1) return 1;
    printf("Limite superior b: ");
    if (scanf("%lf", &b) != 1) return 1;
    if (b <= a) {
        printf("Erro: e necessario que b > a.\n");
        return 1;
    }
    printf("Numero de subintervalos n: ");
    if (scanf("%d", &n) != 1 || n < 1) {
        printf("Erro: n deve ser um inteiro >= 1.\n");
        return 1;
    }

    h = (b - a) / n;
    M = max_f2(a, b);

    if (n == 1) {
        printf("\n>> n = 1: usando a formula de TRAPEZIOS SIMPLES\n");
        I = trapezio_simples(a, b);
        erro = pow(b - a, 3) / 12.0 * M;
    } else {
        printf("\n>> n = %d: usando a formula de TRAPEZIOS GENERALIZADA\n", n);
        I = trapezio_generalizado(a, b, n);
        erro = (b - a) * h * h / 12.0 * M;
    }

    printf("h                    = %.10f\n", h);
    printf("Integral aproximada  = %.10f\n", I);
    printf("max|f''(x)| em [a,b] = %.10f\n", M);
    printf("Margem de erro (|E|<=) = %.10e\n", erro);
    printf("Intervalo estimado   : [%.10f , %.10f]\n", I - erro, I + erro);

#if CONHECE_PRIMITIVA
    {
        double exata = F(b) - F(a);
        printf("\nValor exato          = %.10f\n", exata);
        printf("Erro real absoluto   = %.10e\n", fabs(exata - I));
    }
#endif

    return 0;
}