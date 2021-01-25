#include <stdio.h>

void main()
{
    double Tf;
    printf("Enter temperature in Fahrenheit: ");
    scanf("%lf", &Tf);

    double Tc = (Tf-32.)*(5./9.);
    printf("%.1lf degrees F corresponds to %.1lf degrees C\n", Tf, Tc);
}