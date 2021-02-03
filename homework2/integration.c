#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double func(double x) // to make it easy to switch out the function we're integrating
{
    return exp(x);
}

double integrate_bins(double a, double b, int Nbins) // numerical integrator if you care about setting # of bins
{
    double dx = (b-a)/Nbins;
    double sum = 0;
    for(int i=0; i<Nbins; i++)
    {
        sum += func(a+i*dx+0.5*dx)*dx;
    }
    return sum;
}

double integrate_steps(double a, double b, double dx) // numerical integrator if you care about a specific bin width
{
    double sum = 0;
    for (double x=a; x<b; x+=dx)
    {
        sum += func(x+(0.5*dx))*dx;
    }
    return sum;
}

int main(int argc, char *argv[]) // pass # of bins argument from command line
{
    int bins = atoi(argv[1]);
    double integral = integrate_bins(0,1,bins);
    printf("\nThe integral from 0 to 1 of e^x using %d bins is %lf\n",bins,integral);

    return 0;
}