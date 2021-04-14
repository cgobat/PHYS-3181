#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "integrator.h"

int N;
#define DIM 2

// compute force based on Lennard-Jones potential
void compute_force(double *ri, double *rj, double *fij)
{
    double rij[DIM];
    for (int a = 0; a < DIM; ++a)
        rij[a] = rj[a] - ri[a];
    double rsq = 0;
    for (int a = 0; a < DIM; ++a)
        rsq += rij[a] * rij[a];
    double r = sqrt(rsq);
    for (int a = 0; a < DIM; ++a)
        fij[a] = -24 * (2 / pow(r, 12) - 1 / pow(r, 6)) * rij[a] / (r * r);
}

// the layout for the array x
// (x0 y0 vx0 vy0, x1 y1 vx1 vy1, ...)
void deriv(double *x, double t, double *res)
{
    for (int i = 0; i < N; ++i)
    {
        for (int a = 0; a < DIM; ++a)
            res[i * 2 * DIM + a] = x[i * 2 * DIM + a + DIM];

        double totalforce[DIM];
        for (int a = 0; a < DIM; ++a)
            totalforce[a] = 0;

        for (int j = 0; j < N; ++j)
            if (i != j)
            {
                double force[DIM];
                compute_force(&x[i * 2 * DIM], &x[j * 2 * DIM], force);
                for (int a = 0; a < DIM; ++a)
                    totalforce[a] += force[a];
            }

        for (int a = 0; a < DIM; ++a)
            res[i * 2 * DIM + a + DIM] = totalforce[a];
    }
}

double lj_potential(double *x, double *y)
{
    double dx[2] = {x[0] - y[0], x[1] - y[1]};
    double r = sqrt(dx[0] * dx[0] + dx[1] * dx[1]);

    return 4 * (pow(r, -12) - pow(r, -6));
}

double potential_energy(double *config, int N)
{
    double res = 0;
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            res += lj_potential(&config[i * 4], &config[j * 4]);

    return res;
}

void initial_state(double *x, double L, double e)
{
    for (int i = 0; i < 2 * N * DIM; ++i)
        x[i] = 0;

    // first set positions on a grid
    int sz = ceil(sqrt(N));
    double a = pow(2.0, 1. / 6);

    if (sz * a > L)
        a = L / sz;

    int n = 0;
    for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j, ++n)
            if (n < N)
            {
                x[4 * n + 0] = i * a;
                x[4 * n + 1] = j * a;
            }

    double poten = potential_energy(x, N);

    double ken = N * e - poten;
    if (ken < 0)
    {
        printf("Energy per particle is too low ...\n");
        abort();
    }

    double v = sqrt(2 * ken / N);
    for (int i = 0; i < N; ++i)
    {
        double ang = 2 * M_PI * drand48();
        x[4 * i + 2] = v * cos(ang);
        x[4 * i + 3] = v * sin(ang);
    }
}

int main(int argc, char **argv)
{
    double dt = atof(argv[1]);
    int nstep = atoi(argv[2]);
    N = atoi(argv[3]);
    double L = atof(argv[4]);
    double e = atof(argv[5]);

    srand48(123);

    //printf("argc: %d dt: %f nstep: %d\n", argc, dt, nstep);

    double *x = malloc(2 * N * DIM * sizeof(double));
    double *xf = malloc(2 * N * DIM * sizeof(double));

    initial_state(x, L, e);

    printf("%e ", 0.0);
    for (int i = 0; i < 2 * N * DIM; ++i)
        printf("%e ", x[i]);
    printf("\n");

    for (int i = 0; i < nstep; ++i)
    {
        rk2(x, 2 * N * DIM, i * dt, dt, deriv, xf);

        // check whether particles are still in the box
        for (int k = 0; k < N; ++k)
        {
            if (xf[4 * k] < 0)
            {
                xf[4 * k + 0] *= -1;
                xf[4 * k + 2] *= -1;
            }
            if (xf[4 * k] > L)
            {
                xf[4 * k + 0] = 2 * L - xf[4 * k];
                xf[4 * k + 2] *= -1;
            }
            if (xf[4 * k + 1] < 0)
            {
                xf[4 * k + 1] *= -1;
                xf[4 * k + 3] *= -1;
            }
            if (xf[4 * k + 1] > L)
            {
                xf[4 * k + 1] = 2 * L - xf[4 * k + 1];
                xf[4 * k + 3] *= -1;
            }
        }

        printf("%e ", i * dt + dt);
        for (int j = 0; j < 2 * N * DIM; ++j)
            printf("%e ", xf[j]);
        printf("\n");
        for (int j = 0; j < 2 * N * DIM; ++j)
            x[j] = xf[j];
    }

    return 0;
}
