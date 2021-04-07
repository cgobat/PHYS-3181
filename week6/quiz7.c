#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void compute_stats(double* x, int N, double* av, double* std)
{
    double result = 0;
    double res_sq = 0;
    for(int i=0; i<N; i++) {result += x[i]; res_sq += x[i]*x[i];}
    *av = result/N;
    *std = sqrt(res_sq/N - pow(*av,2));
}

double gaussrnd()
{
    double u1 = drand48();
    double u2 = drand48();
    return -sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

int main()
{
    double sample[100];
    for(int i=0; i<100; i++) sample[i] = 1 + 2*gaussrnd();

    double mean, stdev;
    compute_stats(sample,100,&mean,&stdev);
    printf("Mean: %lf\nStdev: %lf\n",mean,stdev);
    return 0;
}