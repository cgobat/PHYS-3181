#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalgebra.h"

struct dpoint{double x,y, sigma;};
struct trial{int n; double *x; double *y;};

double gaussrnd()
{
  double u1 = drand48();
  double u2 = drand48();
  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

void generate_data(struct dpoint* data, int N)
{
    for(int i=0; i<N; i++)
    {
        data[i].x = i;
        data[i].y = 2*i + i*i + 3 * drand48()-0.5;
        data[i].sigma = 0.5;
    }
}

double poly(double* alpha, int M, double x) // polynomial of degree M
{
    double res = 0;
    for(int k=0; k<M+1; k++) res += alpha[k]*pow(x, k);
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

void stats(double* y, int N, double *av, double *std) // statistics for a single array of datapoints
{
  double sum=0, sumsq=0;
    for(int i=0; i<N; i++)
    {
      sum += y[i];
      sumsq += y[i]*y[i];
    }
  *av = sum/N;
  sumsq /= N;

  *std = sqrt(sumsq-(*av)*(*av));
}

void compute_stats(struct trial* y, int N, int n, double* av, double* std) // like above, but for an array of arrays
{
  double sum=0, sumsq=0;
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<N; j++)
        {
          sum += y[i].y[j];
          sumsq += y[i].y[j]*y[i].y[j];
        }
      av[i] = sum/N;
      sumsq /= N;
      std[i] = sqrt(sumsq-(av[i])*(av[i]));
      sum=0;
      sumsq=0;
      }
}

void resample(struct dpoint* original, struct dpoint* new, int N) // resample value of the y-coordinate within Gaussian error
{
   for(int i=0; i<N; i++)
   {
      new[i].y=original[i].y + original[i].sigma*gaussrnd();
      new[i].x=original[i].x;
      new[i].sigma=original[i].sigma;
   }
}

double coefficients(struct trial* entries, int N_entries, int i) // compute the polynomial coefficient at index i
{
    double *res = malloc(N_entries*sizeof(double));
    int N_points = 13;
    struct dpoint* all = malloc(N_points*sizeof(struct dpoint));

    for(int i=0; i<N_points; i++)
    {
        for(int j=0; j<N_entries; j++)
        {
            res[j] = entries[j].y[i];
            stats(res, N_entries, &all[i].y, &all[i].sigma);
            all[i].x = entries[0].x[i];
        }
    }

    double a[3]; // quadratic -> 3 terms
    chisq_fit(all, N_points, a, 3);
    free(all);
    free(res);
    return a[i];
}

void bootstrap(struct trial* y, int N, int k, double *average, double* error, int coef_i) // perform bootstrap routine
{
    double *result = malloc(k*sizeof(double));
    struct trial* array = malloc(N*sizeof(struct trial));

    for(int i=0; i<k; i++)
    {
        for(int j=0; j<N; j++) array[j] = y[(int)(floor)(drand48()*N)]; // get value at a random integer index
        result[i] = coefficients(array, N, coef_i);
    }

    stats(result, k, average, error);

    free(array);
    free(result);
}

int main(int argc, char *argv[])
{

    int num_resamples = atoi(argv[2]);

    FILE *filepointer;

    if((filepointer = fopen(argv[1], "r")) ==  NULL){
    printf("Failed to open file...");
    exit(1);
    }

    int N, n;
    fscanf(filepointer, "%d %d", &N, &n); // first line
    struct trial* y = malloc(N*sizeof(struct trial)); // array of trial objects, length N
    struct trial* b = malloc(n*sizeof(struct trial)); // array of trial objects, length n
    struct dpoint* c = malloc(n*sizeof(struct dpoint)); // array of datapoints, length n
    struct dpoint* other = malloc(n*sizeof(struct dpoint)); // array of datapoints, length n

    // initialize
    y[0].n = n;
    b[0].n = n;

    for(int i=0; i<N; i++)
    {
        y[i].x = malloc(n*sizeof(double));
        y[i].y = malloc(n*sizeof(double));
    }

    for(int i=0; i<n; i++)
    {
        b[i].x = malloc(N*sizeof(double));
        b[i].y = malloc(N*sizeof(double));
    }

    for(int i=0; i<n; i++)
    { // now read the rest of the file
        for(int j=0; j<N; j++) fscanf(filepointer, "%lf %lf", &y[j].x[i], &y[j].y[i]);
    }

    fclose(filepointer);

    for(int i = 0; i<n; i++)
    {
        for(int j=0; j<N; j++)
        { // correlate y-values
            b[i].y[j] = y[j].y[i];
        }
    }

    double* mean = malloc(n*sizeof(double));
    double* stdev = malloc(n*sizeof(double));

    compute_stats(b, N, n, mean, stdev);

    for(int i=0; i<n; i++)
    {
        c[i].y=mean[i];
        c[i].x=.5*(i+1);
        c[i].sigma=stdev[i];
    }

    for(int i=0;i<num_resamples;i++) resample(c, other, n);
    //for(int i=0; i<n; i++) printf("%e %e\n", other[i].x, other[i].y); //output file contents

    for(int i=0; i<3; i++)
    {
      bootstrap(y, N, n, mean, stdev, i);
      printf("Coefficient %d bootstrapped mean: %e (std dev: %e)\n", i+1, mean, stdev);
    }

    free(y);
    free(b);
    free(c);
    free(other);
    free(mean);
    free(stdev);
    return 0;
}
