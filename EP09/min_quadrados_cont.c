/*
 * CMAC05 - Calculo Numerico para Computacao
 * EP09 - Metodo dos Minimos Quadrados (MMQ)
 *
 * Programa generico (C99 padrao, sem bibliotecas externas alem da libm).
 *
 * Modo 1: MMQ CONTINUO - aproxima uma funcao f(x) em [a,b] por um
 *         polinomio de grau n (qualquer funcao, intervalo e grau).
 *         Basta editar a funcao f() logo abaixo.
 *         Sistema normal:  sum_j c_j <x^i, x^j> = <f, x^i>
 *         <x^i,x^j> = (b^(i+j+1) - a^(i+j+1)) / (i+j+1)   (exato)
 *         <f,x^i>   = integral de f(x)*x^i em [a,b]       (Simpson)
 *
 * Modo 2: MMQ DISCRETO com modelo logistico  g(t) = L / (1 + alpha*exp(-beta*t))
 *         (L informado pelo usuario, pontos (t,F) quaisquer).
 *         Linearizacao: ln(L/F - 1) = ln(alpha) - beta*t  -> reta por MMQ.
 *
 * Compilar:  gcc -std=c99 -Wall -O2 ep09_mmq.c -o ep09_mmq -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXG   12      /* grau maximo + 1 (matrizes ate 12x12) */
#define NSIMP  4000    /* subintervalos (par) da regra de Simpson */

/* ------------------------------------------------------------------ */
/* >>> EDITE AQUI a funcao do Modo 1 <<<                              */
/* Exercicio 1: f(x) = (x+1)/(x^2+1)                                  */
/* ------------------------------------------------------------------ */
static double f(double x)
{
    return (x + 1.0) / (x * x + 1.0);
}

/* ------------------------------------------------------------------ */
/* Utilitarios numericos                                              */
/* ------------------------------------------------------------------ */

/* Resolve A x = b por eliminacao de Gauss com pivotamento parcial.
 * Retorna 0 se OK, -1 se a matriz for (numericamente) singular. */
static int resolve_sistema(int n, double A[MAXG][MAXG], double b[MAXG],
                           double x[MAXG])
{
    int i, j, k, p;
    double tmp, m;

    for (k = 0; k < n - 1; k++) {
        p = k;
        for (i = k + 1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[p][k])) p = i;
        if (fabs(A[p][k]) < 1e-14) return -1;
        if (p != k) {
            for (j = 0; j < n; j++) {
                tmp = A[k][j]; A[k][j] = A[p][j]; A[p][j] = tmp;
            }
            tmp = b[k]; b[k] = b[p]; b[p] = tmp;
        }
        for (i = k + 1; i < n; i++) {
            m = A[i][k] / A[k][k];
            for (j = k; j < n; j++) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }
    if (fabs(A[n - 1][n - 1]) < 1e-14) return -1;

    for (i = n - 1; i >= 0; i--) {
        tmp = b[i];
        for (j = i + 1; j < n; j++) tmp -= A[i][j] * x[j];
        x[i] = tmp / A[i][i];
    }
    return 0;
}

/* Avalia polinomio c[0] + c[1]x + ... + c[n]x^n (Horner) */
static double poly_eval(const double c[MAXG], int n, double x)
{
    double r = c[n];
    int i;
    for (i = n - 1; i >= 0; i--) r = r * x + c[i];
    return r;
}

/* Integrando de <f, x^i> */
static double integrando(double x, int i)
{
    return f(x) * pow(x, i);
}

/* Simpson composta de f(x)*x^i em [a,b] */
static double simpson_fxi(double a, double b, int i)
{
    double h = (b - a) / NSIMP, s = integrando(a, i) + integrando(b, i);
    int k;
    for (k = 1; k < NSIMP; k++)
        s += (k % 2 ? 4.0 : 2.0) * integrando(a + k * h, i);
    return s * h / 3.0;
}

/* Simpson composta de (f - p)^2 em [a,b] (norma L2 do erro) */
static double erro_quadratico(double a, double b, const double c[MAXG], int n)
{
    double h = (b - a) / NSIMP, s, e;
    int k;
    e = f(a) - poly_eval(c, n, a); s = e * e;
    e = f(b) - poly_eval(c, n, b); s += e * e;
    for (k = 1; k < NSIMP; k++) {
        double x = a + k * h;
        e = f(x) - poly_eval(c, n, x);
        s += (k % 2 ? 4.0 : 2.0) * e * e;
    }
    return s * h / 3.0;
}

/* MMQ discreto: ajusta polinomio de grau deg a m pontos (x,y).
 * Retorna 0 se OK. */
static int polyfit_discreto(const double *x, const double *y, int m, int deg,
                            double c[MAXG])
{
    double A[MAXG][MAXG], b[MAXG], pot[2 * MAXG];
    int i, j, k;

    for (i = 0; i <= 2 * deg; i++) {
        pot[i] = 0.0;
        for (k = 0; k < m; k++) pot[i] += pow(x[k], i);
    }
    for (i = 0; i <= deg; i++) {
        for (j = 0; j <= deg; j++) A[i][j] = pot[i + j];
        b[i] = 0.0;
        for (k = 0; k < m; k++) b[i] += y[k] * pow(x[k], i);
    }
    return resolve_sistema(deg + 1, A, b, c);
}

/* ------------------------------------------------------------------ */
/* Modo 1: MMQ continuo polinomial                                    */
/* ------------------------------------------------------------------ */
static void modo_continuo(void)
{
    double a, b, A[MAXG][MAXG], rhs[MAXG], c[MAXG];
    int n, i, j, m, k;
    FILE *fd, *fg;

    printf("Funcao f(x) definida no codigo (edite f() para trocar).\n");
    printf("Intervalo [a b]: ");
    if (scanf("%lf %lf", &a, &b) != 2 || !(a < b)) {
        printf("Intervalo invalido.\n"); return;
    }
    printf("Grau do polinomio (0 a %d): ", MAXG - 1);
    if (scanf("%d", &n) != 1 || n < 0 || n >= MAXG) {
        printf("Grau invalido.\n"); return;
    }

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++)
            A[i][j] = (pow(b, i + j + 1) - pow(a, i + j + 1)) / (i + j + 1);
        rhs[i] = simpson_fxi(a, b, i);
    }

    printf("\nSistema normal (matriz | termo independente):\n");
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) printf("%12.6f ", A[i][j]);
        printf("| %12.6f\n", rhs[i]);
    }

    if (resolve_sistema(n + 1, A, rhs, c) != 0) {
        printf("Sistema singular.\n"); return;
    }

    printf("\nPolinomio h%d(x) = ", n);
    for (i = 0; i <= n; i++) printf("%+.8f*x^%d ", c[i], i);
    printf("\n");
    printf("Norma L2 do erro ||f-h||_2 = %.8f\n",
           sqrt(erro_quadratico(a, b, c, n)));

    printf("\nQuantos pontos deseja avaliar? ");
    if (scanf("%d", &m) != 1 || m < 0) return;
    printf("\n%12s %14s %14s %14s\n", "x", "f(x)", "h(x)", "|f-h|");
    for (k = 0; k < m; k++) {
        double x, fx, hx;
        printf("x[%d] = ", k + 1);
        if (scanf("%lf", &x) != 1) return;
        fx = f(x); hx = poly_eval(c, n, x);
        printf("%12.6f %14.8f %14.8f %14.8f\n", x, fx, hx, fabs(fx - hx));
    }

    /* Dados para grafico (gnuplot) */
    fd = fopen("grafico_modo1.dat", "w");
    if (fd) {
        for (k = 0; k <= 400; k++) {
            double x = a + (b - a) * k / 400.0;
            fprintf(fd, "%.8f %.8f %.8f\n", x, f(x), poly_eval(c, n, x));
        }
        fclose(fd);
        fg = fopen("grafico_modo1.gp", "w");
        if (fg) {
            fprintf(fg, "set grid\nset xlabel 'x'\nset ylabel 'y'\n"
                        "plot 'grafico_modo1.dat' u 1:2 w l lw 2 t 'f(x)', "
                        "'' u 1:3 w l lw 2 dt 2 t 'h(x)'\npause -1\n");
            fclose(fg);
        }
        printf("\nArquivos gerados: grafico_modo1.dat e grafico_modo1.gp\n"
               "(execute:  gnuplot grafico_modo1.gp)\n");
    }
}

/* ------------------------------------------------------------------ */
/* Modo 2: MMQ discreto, modelo logistico g(t)=L/(1+alpha e^(-beta t)) */
/* ------------------------------------------------------------------ */
static void modo_logistico(void)
{
    double L, *t, *F, *y, c[MAXG], alpha, beta, sse = 0.0;
    int m, k;
    FILE *fd, *fg;

    printf("Valor de saturacao L (ex.: 100): ");
    if (scanf("%lf", &L) != 1 || L <= 0.0) { printf("L invalido.\n"); return; }
    printf("Numero de pontos: ");
    if (scanf("%d", &m) != 1 || m < 2) { printf("Pontos insuficientes.\n"); return; }

    t = malloc(m * sizeof *t); F = malloc(m * sizeof *F); y = malloc(m * sizeof *y);
    if (!t || !F || !y) { printf("Sem memoria.\n"); free(t); free(F); free(y); return; }

    for (k = 0; k < m; k++) {
        printf("t[%d] F[%d] = ", k + 1, k + 1);
        if (scanf("%lf %lf", &t[k], &F[k]) != 2) goto fim;
        if (F[k] <= 0.0 || F[k] >= L) {
            printf("F deve estar em (0, L) para linearizar.\n"); goto fim;
        }
        y[k] = log(L / F[k] - 1.0);    /* y = ln(alpha) - beta*t */
    }

    if (polyfit_discreto(t, y, m, 1, c) != 0) {
        printf("Sistema singular.\n"); goto fim;
    }
    alpha = exp(c[0]);
    beta  = -c[1];

    printf("\nLinearizacao: y = ln(L/F - 1) = a0 + a1*t\n");
    printf("a0 = %.8f   a1 = %.8f\n", c[0], c[1]);
    printf("alpha = exp(a0) = %.8f\nbeta  = -a1     = %.8f\n", alpha, beta);
    printf("\ng(t) = %.4f / (1 + %.6f * exp(-%.6f * t))\n", L, alpha, beta);

    printf("\n%10s %12s %12s %12s\n", "t", "F(t)", "g(t)", "F-g");
    for (k = 0; k < m; k++) {
        double g = L / (1.0 + alpha * exp(-beta * t[k]));
        printf("%10.4f %12.4f %12.4f %12.4f\n", t[k], F[k], g, F[k] - g);
        sse += (F[k] - g) * (F[k] - g);
    }
    printf("Soma dos quadrados dos residuos (escala original) = %.6f\n", sse);

    fd = fopen("grafico_modo2.dat", "w");
    if (fd) {
        double tmin = t[0], tmax = t[0];
        for (k = 1; k < m; k++) {
            if (t[k] < tmin) tmin = t[k];
            if (t[k] > tmax) tmax = t[k];
        }
        for (k = 0; k <= 400; k++) {
            double tt = tmin + (tmax - tmin) * k / 400.0;
            fprintf(fd, "%.8f %.8f\n", tt, L / (1.0 + alpha * exp(-beta * tt)));
        }
        fclose(fd);
        fd = fopen("pontos_modo2.dat", "w");
        if (fd) {
            for (k = 0; k < m; k++) fprintf(fd, "%.8f %.8f\n", t[k], F[k]);
            fclose(fd);
        }
        fg = fopen("grafico_modo2.gp", "w");
        if (fg) {
            fprintf(fg, "set grid\nset xlabel 't'\nset ylabel 'F'\n"
                        "plot 'pontos_modo2.dat' w p pt 7 t 'dados', "
                        "'grafico_modo2.dat' w l lw 2 t 'g(t)'\npause -1\n");
            fclose(fg);
        }
        printf("\nArquivos gerados: grafico_modo2.dat, pontos_modo2.dat, "
               "grafico_modo2.gp\n");
    }
fim:
    free(t); free(F); free(y);
}

/* ------------------------------------------------------------------ */
int main(void)
{
    int op;
    printf("=== EP09 - Metodo dos Minimos Quadrados ===\n");
    printf("1) MMQ continuo: aproximar f(x) por polinomio em [a,b]\n");
    printf("2) MMQ discreto: ajuste g(t) = L/(1+alpha*exp(-beta*t))\n");
    printf("Opcao: ");
    if (scanf("%d", &op) != 1) return 1;
    if (op == 1)      modo_continuo();
    else if (op == 2) modo_logistico();
    else              printf("Opcao invalida.\n");
    return 0;
}