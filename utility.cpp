#include "Header.h"

double f(double x, double y, int k)
{
    if (k == 0) 
    {
        return 1;
    }
    else if (k == 1)
    {
        return x;
    }
    else if (k == 2)
    {
        return y;
    }
    else if (k == 3)
    {
        return x + y;
    }
    else if (k == 4)
    {
        return sqrt(x*x + y*y);
    }
    else if (k == 5)
    {
        return x*x + y*y;
    }
    else if (k == 6)
    {
        return exp(x*x - y*y);
    }
    else if (k == 7)
    {
        return 1/(25*(x*x + y*y) + 1);
    }
    return 0;
}

double fxx(double x, double y, int k)
{
    if (k == 0) return 0;
    if (k == 1) return 0;
    if (k == 2) return 0;
    if (k == 3) return 0;

    if (k == 4)
    {
        double r = sqrt(x*x + y*y);
        if (r == 0) return 0;
        return y*y / (r*r*r);
    }

    if (k == 5) return 2;

    if (k == 6)
    {
        return (2 + 4*x*x) * exp(x*x - y*y);
    }

    if (k == 7)
    {
        double q = 25*(x*x + y*y) + 1;
        return -50/(q*q) + 5000*x*x/(q*q*q);
    }

    return 0;
}

double fyy(double x, double y, int k)
{
    if (k == 0) return 0;
    if (k == 1) return 0;
    if (k == 2) return 0;
    if (k == 3) return 0;

    if (k == 4)
    {
        double r = sqrt(x*x + y*y);
        if (r == 0) return 0;
        return x*x / (r*r*r);
    }

    if (k == 5) return 2;

    if (k == 6)
    {
        return (-2 + 4*y*y) * exp(x*x - y*y);
    }

    if (k == 7)
    {
        double q = 25*(x*x + y*y) + 1;
        return -50/(q*q) + 5000*y*y/(q*q*q);
    }

    return 0;
}



void solve_tridiagonal(int n, double *a, double *b, double *c, double *r, double *d )
{
    double *alpha = new double[n];
    double *beta = new double[n];

    alpha[0] = -c[0] / b[0];
    beta[0] = r[0] / b[0];

    for (int i = 1; i < n; i++)
    {
        double denom = b[i] + a[i] * alpha[i - 1];

        if (i != n - 1)
            alpha[i] = -c[i] / denom;

        beta[i] = (r[i] - a[i] * beta[i - 1]) / denom;
    }

    d[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        d[i] = alpha[i] * d[i + 1] + beta[i];
    }

    delete[] alpha;
    delete[] beta;
}

double spline_1d(double *x, double *val, int n, double f2_left, double f2_right, double xp)
{
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *r = new double[n];
    double *d = new double[n];

    // Левая граница:
    // 2d1 + d2 = 3f[x1;x2] - 1/2 f''(x1)(x2-x1)
    double h = x[1] - x[0];
    double div = (val[1] - val[0]) / h;

    a[0] = 0;
    b[0] = 2;
    c[0] = 1;
    r[0] = 3 * div - 0.5 * f2_left * h;

    // Внутренние уравнения
    for (int i = 1; i <= n - 2; i++)
    {
        double h_left = x[i] - x[i - 1];
        double h_right = x[i + 1] - x[i];

        double div_left = (val[i] - val[i - 1]) / h_left;
        double div_right = (val[i + 1] - val[i]) / h_right;

        a[i] = h_right;
        b[i] = 2 * (x[i + 1] - x[i - 1]);
        c[i] = h_left;

        r[i] = 3 * div_left * h_right + 3 * div_right * h_left;
    }

    // Правая граница:
    // d_{n-1} + 2d_n = 3f[x_{n-1};x_n] + 1/2 f''(x_n)(x_n-x_{n-1})
    h = x[n - 1] - x[n - 2];
    div = (val[n - 1] - val[n - 2]) / h;

    a[n - 1] = 1;
    b[n - 1] = 2;
    c[n - 1] = 0;
    r[n - 1] = 3 * div + 0.5 * f2_right * h;

    solve_tridiagonal(n, a, b, c, r, d);

    // Находим отрезок [x[i], x[i+1]]
    int idx = 0;

    if (xp <= x[0])
        idx = 0;
    else if (xp >= x[n - 1])
        idx = n - 2;
    else
    {
        for (int i = 0; i < n - 1; i++)
        {
            if (xp >= x[i] && xp <= x[i + 1])
            {
                idx = i;
                break;
            }
        }
    }

    h = x[idx + 1] - x[idx];
    div = (val[idx + 1] - val[idx]) / h;

    double c1 = val[idx];
    double c2 = d[idx];
    double c3 = (3 * div - 2 * d[idx] - d[idx + 1]) / h;
    double c4 = (d[idx] + d[idx + 1] - 2 * div) / (h * h);

    double t = xp - x[idx];

    double result = c1 + c2 * t + c3 * t * t + c4 * t * t * t;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] r;
    delete[] d;

    return result;
}

