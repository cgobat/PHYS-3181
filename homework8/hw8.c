#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "integrator.h"
#define DIM 2

int N;

// compute force based on Lennard-Jones potential
void compute_force(double* ri, double *rj, double* fij)
{
    double rij[DIM];
    for(int a=0; a<DIM; ++a) rij[a] = rj[a]-ri[a];
    double rsq=0;
    for(int a=0; a<DIM; ++a) rsq += rij[a]*rij[a];
    double r = sqrt(rsq);
    for(int a=0; a<DIM; ++a) fij[a] = -24*(2/pow(r,12)-1/pow(r,6))*rij[a]/(r*r);
}

// the layout for the array x
// (x0 y0 vx0 vy0, x1 y1 vx1 vy1, ...)
void deriv(double *x, double t, double* res)
{
    for(int i=0; i<N; i++)
    {
        for(int a=0; a<DIM; ++a) res[i*2*DIM+a] = x[i*2*DIM+a+DIM];

        double totalforce[DIM];
        for(int a=0; a<DIM; ++a) totalforce[a] = 0;
     
        for(int j=0; j<N; j++) if(i!=j)
        {
            double force[DIM];
            compute_force(&x[i*2*DIM], &x[j*2*DIM], force);
            for(int a=0; a<DIM; ++a) totalforce[a] += force[a];
        }

        for(int a=0; a<DIM; ++a) res[i*2*DIM+a+DIM] = totalforce[a];

    }
}

double lj_potential(double* x, double* y)
{
    double dx[2] = {x[0]-y[0], x[1]-y[1]};
    double r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]);

    return 4*(pow(r,-12)-pow(r,-6));
}
 
double    potential_energy(double* config, int N)
{
    double res = 0;
    for(int i=0; i<N; i++)    
    for(int j=i+1; j<N; j++)
        res += lj_potential(&config[i*4], &config[j*4]);
 
    return res;
}

void initial_state(double *x, double L, double e)
{
    for(int i=0; i<2*N*DIM; i++) x[i] = 0;


    // first set positions on a grid
    int sz = ceil(sqrt(N));
    double a = pow(2.0, 1./6);

    if(sz*a > L) a = L/sz;

    int n=0;
    for(int i=0; i<sz; i++)
    for(int j=0; j<sz; j++, ++n) if(n<N)
    {
        x[4*n+0] = i*a;
        x[4*n+1] = j*a;
    }

    double poten = potential_energy(x, N);

    double ken = N*e - poten;
    if(ken < 0)
    {
        printf("Energy per particle is too low ...\n");
        abort();
    }

    double v = sqrt(2*ken/N);
    for(int i=0; i<N; i++)
    {
        double ang = 2*M_PI*drand48();
        x[4*i+2] = v*cos(ang);
        x[4*i+3] = v*sin(ang);
    }
}

int main(int argc, char *argv[])
{
    N = atoi(argv[1]);
    double L = atof(argv[2]);
    double EperN = atof(argv[3]);
    double dt = atof(argv[4]);
    double t = atof(argv[5]);
    int nskip = atoi(argv[6]);
    int nstep = t/dt;

    printf("%d particle(s), box dimension: %.2lf, "
        "time: %.2lf, dt: %.2lf, %d steps\n",N,L,t,dt,nstep);

    double *x = malloc(2*N*DIM*sizeof(double));
    double *xn = malloc(2*N*DIM*sizeof(double));

    // wall impulses
    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    double a4 = 0;

    initial_state(x, L, EperN);

    printf("%.2lf ", 0);
    for(int i=0; i<2*N*DIM; i++) printf("%.3lf ", x[i]);
    printf("%.3lf %.3lf %.3lf %.3lf\n", a1, a2, a3, a4);

    for(int i=0; i<nstep; i++)
    {
        rk2(x, 2*N*DIM, i*dt, dt, deriv, xn);

        for(int j=0; j<N; j++)
        {
            if(xn[4*j] < 0)
            { 
                xn[4*j+0] *= -1;
                xn[4*j+2] *= -1;
                a1 += 2*xn[4*j+2];    
            }
            if(xn[4*j] > L)
            { 
                xn[4*j+0] = 2*L - xn[4*j];
                xn[4*j+2] *= -1;
                a2 += 2*xn[4*j+2];
            }
            if(xn[4*j+1] < 0) 
            { 
                xn[4*j+1] *= -1; 
                xn[4*j+3] *= -1; 
                a3 += 2*xn[4*j+3];
            }
            if(xn[4*j+1] > L) 
            { 
                xn[4*j+1] = 2*L - xn[4*j+1]; 
                xn[4*j+3] *= -1; 
                a4 += 2*xn[4*j+3];
            }
        }

        if (i%nskip==0)
        {
                printf("%.2lf ", i*dt+dt);
                for(int j=0; j<N; j++)
                    printf("%.3lf %.3lf %.3lf %.3lf ", xn[4*j], xn[4*j+1], xn[4*j+2], xn[4*j+3]);
                printf("%.3lf %.3lf %.3lf %.3lf\n", a1, a2, a3, a4);
        }

        for(int j=0; j<2*N*DIM; j++) x[j] = xn[j];

    }

    return 0;
}
