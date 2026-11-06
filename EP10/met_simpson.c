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

// 2) Quarta derivada f''''(x) (necessária para o cálculo do erro) 
double f4(double x)
{
    return 24;                
}

/* 3) (OPCIONAL) Primitiva F(x), para comparar com o valor exato.
 *    Se não souber a primitiva, defina CONHECE_PRIMITIVA como 0.     */
#define CONHECE_PRIMITIVA 1        /* <-- ALTERAR: 1 = sim, 0 = não */

double F(double x)
{
    return exp(x);                 
}

/* Simpson 1/3 simples (n = 2) */
double simpson_simples(double a, double b)
{
    double h = (b - a) / 2.0;
    return (h / 3.0) * (f(a) + 4.0 * f(a + h) + f(b));
}

/* Simpson 1/3 generalizado (n par, n > 2) */
double simpson_generalizado(double a, double b, int n)
{
    double h = (b - a) / n;
    double soma = f(a) + f(b);
    int i;
    for (i = 1; i < n; i++) {
        double res;
        if(i%2) {
            res = 4.0;
        }
        else {
            res = 2.0;
        }
        soma += res * f(a + i * h);
    }
    return (h / 3.0) * soma;
}

/* Estima max |f''''(x)| em [a,b] por varredura */
double max_f4(double a, double b)
{
    double max = fabs(f4(a)), v, x;
    int i;
    for (i = 1; i <= PONTOS_VARREDURA; i++) {
        x = a + (b - a) * i / PONTOS_VARREDURA;
        v = fabs(f4(x));
        if (v > max){
            max = v;
        }
    }
    return max;
}

int main(void)
{
    double a, b, h, I, M, erro;
    int n;

    printf("=== Integracao pelo Metodo 1/3 de Simpson ===\n");
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

    // Simpson exige n par: se for impar, soma 1
    if (n % 2 != 0) {
        n++;
        printf("n era impar; ajustado para n = %d\n", n);
    }

    h = (b - a) / n;
    M = max_f4(a, b);

    if (n == 2) {
        printf("\n>> n = 2: usando a formula de SIMPSON 1/3 SIMPLES\n");
        I = simpson_simples(a, b);
        erro = pow(b - a, 5) / 2880.0 * M;
    } else {
        printf("\n>> n = %d: usando a formula de SIMPSON 1/3 GENERALIZADA\n", n);
        I = simpson_generalizado(a, b, n);
        erro = (b - a) * pow(h, 4) / 180.0 * M;
    }

    printf("h                      = %.10f\n", h);
    printf("Integral aproximada    = %.10f\n", I);
    printf("max|f''''(x)| em [a,b] = %.10f\n", M);
    printf("Margem de erro (|E|<=) = %.10e\n", erro);
    printf("Intervalo estimado     : [%.10f , %.10f]\n", I - erro, I + erro);

#if CONHECE_PRIMITIVA
    {
        double exata = F(b) - F(a);
        printf("\nValor exato            = %.10f\n", exata);
        printf("Erro real absoluto     = %.10e\n", fabs(exata - I));
    }
#endif

    return 0;
}