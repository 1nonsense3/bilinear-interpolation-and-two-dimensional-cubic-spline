#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <vector>

using namespace std;

#define fxy(i,j) fxy[(i)*ny + (j)]


double f(double x, double y, int k);
void aprox(int n, double *x, double *fx, double *a, double *d);
double lin_ap_val(int nx, int ny, double * x, double * y, double * fxy, double xp, double yp);
void plot_slices_nodes(double ax, double bx, int nx, int ny, double *x, double *y, double *fxy, int k, int nslices);
double spline_2d(double *x, double *y, double *fxy, int nx, int ny, double xp, double yp, int k);
double spline_1d(double *x, double *val, int n, double f2_left, double f2_right, double xp);
void solve_tridiagonal(int n, double *a, double *b, double *c, double *r, double *d );
double fxx(double x, double y, int k);
double fyy(double x, double y, int k);
void plot_slices_nodes_spline(double ax, double bx,int nx, int ny,double *x, double *y,double *fxy,int k,int nslices);
void plot_3d_compare(double ax, double bx, double ay, double by, int nx, int ny, double *x, double *y, double *fxy, int k);
void plot_3d_linear(
    double ax, double bx,
    double ay, double by,
    int nx, int ny,
    double *x, double *y,
    double *fxy,
    int k
);

void plot_3d_linear_error(
    double ax, double bx,
    double ay, double by,
    int nx, int ny,
    double *x, double *y,
    double *fxy,
    int k
);

void plot_3d_spline_error(
    double ax, double bx,
    double ay, double by,
    int nx, int ny,
    double *x, double *y,
    double *fxy,
    int k
);