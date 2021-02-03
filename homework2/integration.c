#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double func(double x) // to make it easy to switch out the function we're integrating
{
    return exp(x);
}

double integrate_p1(double a, double b, double dx) // numerical integrator if you care about a specific bin width
{
    double sum = 0;
    for (double x=a; x<b; x+=dx)
    {
        sum += func(x)*dx;
    }
    return sum;
}

double integrate_p2(double a, double b, int Nbins) // numerical integrator if you care about setting # of bins
{
    double dx = (b-a)/Nbins;
    double sum = 0;
    for(int i=0; i<Nbins; i++)
    {
        sum += func(a+i*dx)*dx;
    }
    return sum;
}

double integrate_p3(double a, double b, int Nbins) // numerical integrator if you care about setting # of bins
{
    double dx = (b-a)/Nbins;
    double sum = 0;
    for(int i=0; i<Nbins; i++)
    {
        sum += func(a+i*dx+0.5*dx)*dx;
    }
    return sum;
}

double integrate_p4(double theta_0, int Nbins) // numerical integrator if you care about setting # of bins
{
    double rads = M_PI*theta_0/180.;
    double dtheta = 2*rads/Nbins;
    double theta;
    double sum = 0;
    //printf("Expected value for small angle: %lf\n",2*M_PI*sqrt(9.81/9.81));
    //printf("theta_0 is %lf radians\n",rads);
    for(int i=0; i<Nbins; i++)
    {
        theta = (0-rads)+i*dtheta;
        sum += 1/sqrt(2*9.81*(cos(theta+0.5*dtheta)-cos(rads))/9.81)*dtheta;
    }
    return sum;
}

int main(int argc, char *argv[]) // pass arguments from command line
{
    //double N[] = {1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9};
    //double tta_0 = atof(argv[1]);
    // double err = fabs(result-(exp(1)-1))/(exp(1)-1);
    for (double tta_0=5; tta_0<51; tta_0++)
    {
        double result = 2*integrate_p4(tta_0,1e6);
        double error = fabs(result - 2*M_PI*sqrt(9.81/9.81))/fabs(2*M_PI*sqrt(9.81/9.81));
        printf("%lf\t%.16lf\n",tta_0,error);
    }

    return 0;
}