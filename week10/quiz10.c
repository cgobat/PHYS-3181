#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void compute_laplacean(int nphi, double* phi, double* psi, int* ne)
{
    for(int i=0; i<nphi; i++)
    {
        psi[i] = 0;
        for(int j=0; j<4; j++) psi[i] += phi[ne[i*4+j]];
        psi[i] -= 4*phi[i];
    }
}

void print_gradient(int nphi, double* phi, int* ne, double* pos, double a)
{
    for(int i = 0; i<nphi; i++)
    {
        double g_x = 0.5*(phi[ne[2*i]]-phi[ne[2*i+1]])/a;
        double g_y = 0.5*(phi[ne[2*i+2]]-phi[ne[2*i+3]])/a;

        printf("gradient at (%.3e, %.3e): %e, %e\n", pos[2*i], pos[2*i+1], g_x, g_y);
    }
    
}