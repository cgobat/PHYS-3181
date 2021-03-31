#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalgebra.h"

#define ORDER 2 // fit a second-order polynomial


struct dpoint {double x, y, sigma;};

double poly(double* alpha, int M, double x)
{
    double res = 0;
    for(int k=0; k<M+1; ++k) res += alpha[k]*pow(x, k);
    return res;
}


double chisq_fit(struct dpoint* data, int M, double *a, int N) // chi-squared minimizer
{
    double *A = malloc((N+1)*(N+1)*sizeof(double));
    double *b = malloc((N+1)*(N+1)*sizeof(double));

    for(int k=0; k<N+1; k++)
    {
        b[k]=0;
        for(int i=0; i<M; i++) b[k] += data[i].y*pow(data[i].x, k)/pow(data[i].sigma, 2);
    }

    for(int k=0; k<N+1; k++)
    for(int j=0; j<N+1; j++)
    {
        A[k*(N+1)+j] = 0;
        for(int i=0; i<M; i++) A[k*(N+1)+j] += pow(data[i].x, j+k)/pow(data[i].sigma, 2);
    }

    solve(A, N+1, b, a);

    free(A);
    free(b);

    double chi_sq = 0;
    for(int i=0; i<M; i++) chi_sq += pow(data[i].y-poly(a, N, data[i].x),2)/pow(data[i].sigma, 2);
    return chi_sq;
}

int main(int argc, char *argv[])
{
    char *path = argv[1];
    int N = atoi(argv[2]);
    struct dpoint* data = malloc(N*sizeof(struct dpoint));
    FILE *filepointer;

    if((filepointer = fopen(path, "r")) ==  NULL){
        printf("Could not find file in current directory.");
        exit(1);
    }

    for(int i=0; i<N; i++) fscanf(filepointer, "%lf %lf %lf", &data[i].x, &data[i].y, &data[i].sigma);

    fclose(filepointer);

    // (un)comment to toggle output of the data that is read:
    //for(int i=0; i<N; i++) printf("%.1lf %lf %lf\n", data[i].x, data[i].y, data[i].sigma);
    double coefs[ORDER+1];
    double chisq = chisq_fit(data, N, coefs, ORDER);
    // (un)comment to toggle output of parameter interpretations:
    //printf("Seconds: %.1lf\nH0: %e\nV0: %e\ng: %e\n", N/2.0, coefs[0], coefs[1], coefs[2]*-2);
    printf("Best-fit polynomial: ");
    for(int i=ORDER; i>1; i--) printf("(%e)t^%d + ",coefs[i],i);
    if(ORDER>1) printf("(%e)t + ",coefs[1]);
    printf("%e\n",coefs[0]);
    printf("Chi squared: %lf\n", chisq);
    printf("Reduced chi-sq: %.3lf\n", chisq/(N-1-ORDER));

    free(data);
    return 0;
}
