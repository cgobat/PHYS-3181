#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void init_velocities(double* y, int nbody, double mass, double kT)
{
    double v = sqrt(2*kT/mass);
    for(int i=0; i<nbody; i++)
    {
        double angle = 2*M_PI*drand48();
        y[4*i+2] = v*cos(angle);
        y[4*i+3] = v*sin(angle);
    }
}

double lj_potential(double* x, double* y)
{
    double dx[2] = {x[0]-y[0], x[1]-y[1]};
    double r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]);

    return 4*(pow(r,-12)-pow(r,-6));
}

double compute_pot_one(double* y, int nbody, int j)
{
    double u=0;
    for(int i=0; i<nbody; i++) if(i != j)
    {
        //u += lj_potential(y+4*i, y+4*j);
        u += lj_potential(&y[4*i], &y[4*j]);
    }

    return u;
}

void update(double* y, int n, double kT, double L, double deltax, double deltav)
{
    double dE;

    // first update velocities
    for(int i=0; i<n; i++)
    {
        double vx = y[4*i+2] + deltav*(2*drand48()-1);
        double vy = y[4*i+3] + deltav*(2*drand48()-1);

        dE = 0.5*(vx*vx+vy*vy-pow(y[4*i+2],2)-pow(y[4*i+3],2));
        if(exp(-dE/kT) > drand48()) // accept
        {
            y[4*i+2] = vx;
            y[4*i+3] = vy;
        }
    }

    // now update positions
    for(int i=0; i<n; i++)
    {
        double oldE = compute_pot_one(y, n, i);
        double ox = y[4*i];
        double oy = y[4*i+1];
        y[4*i] += deltax*(2*drand48()-1);
        y[4*i+1] += deltax*(2*drand48()-1);

        // if outside box, reject
        if( (y[4*i]<0) || (y[4*i]>L) || (y[4*i+1]<0) || (y[4*i+1]>L) )
        {
            y[4*i]=ox;
            y[4*i+1]=oy;
            continue;
        }

        dE = compute_pot_one(y, n, i) - oldE;

        if(exp(-dE/kT) < drand48()) // reject
        {
            y[4*i] = ox;
            y[4*i+1] = oy;
        }
        
    }
}

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    double kT = atof(argv[3]);
    double deltax = atof(argv[4]);
    double deltav = atof(argv[5]);
    int n_updates = atoi(argv[6]);
    int n_skip = atoi(argv[7]);

    double *y = malloc(4*N*sizeof(double));

    // initialize the positions
    // first set positions on a grid
    int sz = ceil(sqrt(N));
    double a = pow(2., 1./6.);

    if(sz*a > L)
        a = L/sz;

    int n=0;
    for(int i=0; i<sz; i++)
    for(int j=0; j<sz; j++, n++) if(n<N)
    {
        y[4*n+0] = i*a;
        y[4*n+1] = j*a;
    }

    init_velocities(y, N, 1., kT); // unit 1 mass

    for(int i=0; i<n_updates; i++)
    {
        for(int j=0; j<n_skip; j++) update(y, N, kT, L, deltax, deltav);
        printf("%03d ", i);
        for(int j=0; j<4*N; j++) printf("%e ", y[j]);
        printf("\n\n");
    }

    free(y);

    return 0;
}