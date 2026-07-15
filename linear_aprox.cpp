#include "Header.h"


double lin_ap_val(int nx, int ny, double * x, double * y, double * fxy, double xp, double yp)
{
    for (int i = 0; i < nx - 1; i++) {
        for (int j = 0; j < ny - 1; j++) {

            double xmin = min(x[i], x[i + 1]);
            double xmax = max(x[i], x[i + 1]);
            double ymin = min(y[j], y[j + 1]);
            double ymax = max(y[j], y[j + 1]);

            if (xp >= xmin && xp <= xmax && yp >= ymin && yp <= ymax) {

                double t = (xp - x[i]) / (x[i + 1] - x[i]);
                double u = (yp - y[j]) / (y[j + 1] - y[j]);

                return (1.0 - t) * (1.0 - u) * fxy(i, j) + (1.0 - t) * u * fxy(i, j + 1) + t * (1.0 - u) * fxy(i + 1, j) + t * u * fxy(i + 1, j + 1);
            }
        }
    }

    return 0.0;
}