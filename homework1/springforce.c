#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 4.*atan(1) // simple way to use math library function to define pi

void main(int argc, char *argv[])
{
    if (argc != 5) {
        printf("\nError: wrong number of input values provided.\nPlease pass 4 arguments: A_x, A_y, phi_x, and phi_y, separated by spaces.\n\n");
        exit(EXIT_FAILURE); // user must provide 4 parameters in order for proper execution
    }

    // 0th element of argv is the program name
    double A_x = atof(argv[1]);
    double A_y = atof(argv[2]);
    double phi_x = atof(argv[3]);
    double phi_y = atof(argv[4]);
    double k = 1.0; // set these to be whatever is desired
    double m = 1.0; // for simplicity, I take them to be 1
    double w = sqrt(k/m);
    double wt;

    for(wt=0.; wt<=2.*PI; wt=wt+(2.*PI/100.))
    {
        double t = wt/w; // actual time can be calculated because we set/know omega
        double x = A_x*sin(wt+phi_x);
        double y = A_y*sin(wt+phi_y);
        printf("%9lf\t%9lf\t%9lf\n",t,x,y); // print out time and corresponding x, y position
    }
}