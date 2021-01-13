#include <stdio.h>

void main()
{
    double Tf;
    printf(" deg F\t deg C\n");
    for (Tf=0.;Tf<=100.;Tf=Tf+5)
    {
        double Tc = (Tf-32.)*(5./9.);
        printf("%6.2lf\t%6.2lf\n", Tf, Tc);
    }
}