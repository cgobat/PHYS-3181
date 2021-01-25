#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x)
{
    return 1/sqrt(2*(2-sin(x)));
}

void main(int argc, char *argv[])
{
    double a = atof(argv[1]);
    double b = atof(argv[2]);
    int Nbins = atoi(argv[3]);
    double dx = (b-a)/Nbins;

    double sum = 0;
    for(int i=0; i<Nbins; i++) sum += f(a+i*dx+0.5*dx)*dx;

    printf("The integral of f(x) over the interval [%.2f, %.2f] using %i bins is %lf",a,b,Nbins,sum);
}