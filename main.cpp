#include "Header.h"

int main(int argc, char * argv[])
{
    // входные параметры

    double ap = stod(argv[1]);
    double bp = stod(argv[2]);
    double cp = stod(argv[3]);
    double dp = stod(argv[4]);
    int nx = stoi(argv[5]);
    int ny = stoi(argv[6]);
    int mx = stoi(argv[7]);
    int my = stoi(argv[8]);
    int k = stoi(argv[9]);


    double * x = new double[nx]; // равномерно распределяем точки по оси х
    for (int i = 0; i < nx; i++)
    {
        x[i] = ap + i * (bp - ap)/ (nx - 1);
    }

    double * y = new double[ny]; // равномерно распределяем точки по оси y
    for (int i = 0; i < ny; i++)
    {
        y[i] = cp + i * (dp - cp) / (ny - 1);
    }


    double * fxy = new double[nx*ny]; // считаем значение функции в узлах
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fxy(i, j) = f(x[i], y[j], k);
        }
    }
    plot_3d_compare(ap, bp, cp, dp, nx, ny, x, y, fxy, k);
    plot_3d_spline_error(ap, bp, cp, dp, nx, ny, x, y, fxy, k);

    //plot_3d_linear(ap, bp, cp, dp, nx, ny, x, y, fxy, k);
    //plot_3d_linear_error(ap, bp, cp, dp, nx, ny, x, y, fxy, k);
}