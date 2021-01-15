#include <stdio.h>
#include <math.h>

void main()
{
    double A_x = 2.;
    double A_y = 1.;
    double phi_x = 0.;
    double phi_y = phi_x + M_PI_4;
    double k = 1.0;
    double m = 1.0;
    double w = sqrt(k/m);
    double wt;

    //printf("Enter values for A_x, A_y, phi_x, and phi_y separated by commas: ");
    //scanf("%lf,%lf,%lf,%lf",&A_x,&A_y,&phi_x,&phi_y);
    //printf("%9s\t%9s\t%9s\n","t","x","y");

    for(wt=0.; wt<=2.*M_PI; wt=wt+(2.*M_PI/100.))
    {
        double t = wt/w;
        double x = A_x*sin(wt+phi_x);
        double y = A_y*sin(wt+phi_y);
        printf("%9lf\t%9lf\t%9lf\n",t,x,y);
    }
}