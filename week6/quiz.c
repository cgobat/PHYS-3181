#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct dpoint{double x,y, sigma;};

double poly(double* alpha, int M, double x)
{
    double sum;
    for(int k=0; k<M+1; k++) sum += alpha[k]*pow(x,k);
    return sum;
}

double chisq(double* alpha, int M, struct dpoint* data, int N)
{
    double chi_sq = 0;
    for(int i=0; i<N; i++) chi_sq += pow((data[i].y - poly(alpha, M, data[i].x))/data[i].sigma,2);
    return chi_sq;
}

int main()
{
    struct dpoint data[] = {1,1,1, 2,2,1, 3,4,1, 4,3,1, 5,4.5,0.2};
    double alphas[2] = {0,1}; // constant and linear terms only

    printf("Chi-squared: %f",chisq(alphas,1,data,5));
    
    return 0;
}