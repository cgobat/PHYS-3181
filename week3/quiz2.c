// to compile: gcc -o executable_name quiz2.c

# include <stdio.h>
# include <math.h>

int main()
{
    double x;

    for (int i=0; i<101; i++)
    {
        x = 1. + i*0.01;
        printf("%.2lf\t%lf\n",x,1/sqrt(x));
    }

    return 0;
}