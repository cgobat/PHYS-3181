#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double lj_potential(double* r1, double* r2)
{
    double r = sqrt(pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2));
    return 4*(pow(r,-12) - pow(r,-6));
}

double potential_energy(double* config, int N)
{
    double U = 0;
    for(int i=0; i<N; i++)
    {
        double* ri = &config[4*i];
        for(int j=i+1; j<N; j++)
        {
            double* rj = &config[4*j];
            U += lj_potential(ri,rj);
        }
    }
    return U;
}