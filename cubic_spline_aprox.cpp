#include "Header.h"

double spline_2d(double *x, double *y, double *fxy, int nx, int ny, double xp, double yp, int k)
{
    double *z = new double[ny];

    for (int j = 0; j < ny; j++)
    {
        double *values_x = new double[nx];

        for (int i = 0; i < nx; i++)
        {
            values_x[i] = fxy(i, j);
        }

        double f2_left = fxx(x[0], y[j], k);
        double f2_right = fxx(x[nx - 1], y[j], k);

        z[j] = spline_1d(x, values_x, nx, f2_left, f2_right, xp);

        delete[] values_x;
    }

    double f2_bottom = fyy(xp, y[0], k);
    double f2_top = fyy(xp, y[ny - 1], k);

    double result = spline_1d(y, z, ny, f2_bottom, f2_top, yp);

    delete[] z;

    return result;
}