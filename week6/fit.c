#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct dpoint{double x,y, sigma;};

void generate_data(struct dpoint* data, int N)
{
    for(int i=0; i<N; i++)
    {
        data[i].x = i;
        data[i].y = 2*i + i*i + 3 * drand48()-0.5;
        data[i].sigma = 0.5;
    }
}

void fit_poly(struct dpoint* data, int N, double* coef, int ncoef)
{
    double *A = malloc(ncoef*ncoef*sizeof(double));
    double *b = malloc(ncoef*sizeof(double));

    for(int i=0; i<ncoef; i++)
    {
        for(int j=0; i<ncoef; j++)
        {
            double tmp = 0;
            for(int k=0; k<N; k++) tmp += pow(data[k].x,i+j)/pow(data[k].sigma,2);
            A[i*ncoef+j] = tmp;
        }
    }

    for(int i=0; i<ncoef; i++)
    {
        double tmp=0;
        for(int k=0; k<N; k++) tmp += data[i].y*pow(data[i].x,i)/pow(data[k].sigma, 2);
        b[i] = tmp;
    }

    solve(A, ncoef, b, coef);

    free(A);
    free(b);
}

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

int main(int argc, char *argv[])
{
    int N = 100;
    struct dpoint *data = malloc(N*sizeof(struct dpoint));

    srand48(atoi(argv[1]));

    generate_data(data,N);

    for(int i=0;i<N;i++) printf("%.3f %.3f %.3f\n", data[i].x, data[i].y, data[i].sigma);

    double coef[3];
    fit_poly(data, N, coef, 3);

    printf("coef: "); for(int i=0; i<3; i++) printf("%e ",coef[i]); printf("\n");
    printf("Chi-sq: ",chisq(coef,))

    return 0;
}