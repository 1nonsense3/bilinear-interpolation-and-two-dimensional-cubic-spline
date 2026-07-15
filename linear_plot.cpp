#include "Header.h" 

void plot_3d_linear(
    double ax, double bx,
    double ay, double by,
    int nx, int ny,
    double *x, double *y,
    double *fxy,
    int k
)
{
    FILE *gp = popen("gnuplot -persist", "w");

    if (!gp) {
        std::cerr << "Cannot start gnuplot\n";
        return;
    }

    const int samples_x = 100;
    const int samples_y = 100;

    fprintf(gp, "set term qt\n");

    fprintf(gp, "set hidden3d\n");
    fprintf(gp, "set grid\n");

    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set zlabel 'z'\n");

    fprintf(gp,
        "set xrange [%lf:%lf]\n",
        ax, bx);

    fprintf(gp,
        "set yrange [%lf:%lf]\n",
        ay, by);

    fprintf(gp,
        "splot '-' with lines title 'Linear 2D'\n");

    for (int j = 0; j < samples_y; j++)
    {
        double yp =
            ay + (by - ay) * j / (samples_y - 1);

        for (int i = 0; i < samples_x; i++)
        {
            double xp =
                ax + (bx - ax) * i / (samples_x - 1);

            double zp =
                lin_ap_val(
                    nx, ny,
                    x, y,
                    fxy,
                    xp, yp
                );

            fprintf(gp,
                "%lf %lf %lf\n",
                xp, yp, zp);
        }

        fprintf(gp, "\n");
    }

    fprintf(gp, "e\n");

    fflush(gp);
    pclose(gp);
}

void plot_3d_linear_error(
    double ax, double bx,
    double ay, double by,
    int nx, int ny,
    double *x, double *y,
    double *fxy,
    int k
)
{
    FILE *gp = popen("gnuplot -persist", "w");

    if (!gp) {
        std::cerr << "Cannot start gnuplot\n";
        return;
    }

    const int samples_x = 100;
    const int samples_y = 100;

    fprintf(gp, "set term qt\n");

    fprintf(gp, "set hidden3d\n");
    fprintf(gp, "set grid\n");

    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set zlabel '|f-L|'\n");

    fprintf(gp, "set xrange [%lf:%lf]\n", ax, bx);
    fprintf(gp, "set yrange [%lf:%lf]\n", ay, by);

    fprintf(gp, "splot '-' with lines title 'Linear error'\n");

    for (int j = 0; j < samples_y; j++)
    {
        double yp = ay + (by - ay) * j / (samples_y - 1);

        for (int i = 0; i < samples_x; i++)
        {
            double xp = ax + (bx - ax) * i / (samples_x - 1);

            double exact = f(xp, yp, k);

            double approx =
                lin_ap_val(
                    nx, ny,
                    x, y,
                    fxy,
                    xp, yp
                );

            double err = fabs(exact - approx);

            fprintf(gp, "%lf %lf %lf\n", xp, yp, err);
        }

        fprintf(gp, "\n");
    }

    fprintf(gp, "e\n");

    fflush(gp);
    pclose(gp);
}