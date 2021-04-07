#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void compute_stats(double* x, int N, double *av, double *std)
{
  double sum=0, sumsq=0;
  for(int i=0; i<N; ++i)
  {
    sum += x[i];
    sumsq += x[i]*x[i];
  }

  *av = sum/N;
  *std = sqrt(sumsq/N - (*av)*(*av));

}

double gaussrnd()
{
  double u1 = drand48();
  double u2 = drand48();
  return -sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

void generate_data(double* x, int N)
{
  for(int i=0; i<N; ++i) x[i] = 1 + 2*gaussrnd();
}

void bootstrap(double* x, int N, double* avf, double* erf)
{
  int nboot = 99;
  double *f = malloc(nboot*sizeof(double));

  double *sample = malloc(N*sizeof(double));

  for(int i=0; i<nboot; ++i)
  {
    for(int j=0; j<N; ++j)
    {
      int k = floor(drand48()*N);
      sample[j] = x[k];
    }
    double dummy;
    compute_stats(sample, N, &dummy, &f[i]);
  }

  compute_stats(f, nboot, avf, erf);

  free(sample);
  free(f);
}

int main(int argc, char** argv)
{
  int seed = (argc>1) ? atoi(argv[1]) : 0;
  srand48(seed);

  int N=100;
  double *x = malloc(N*sizeof(double));
  generate_data(x, N);

  double av, std;
  compute_stats(x, N, &av, &std);
  printf("av: %e  std: %e \n", av, std);

  bootstrap(x, N, &av, &std);
  //printf("mean(witherr): %e %e\n", av, std);
  printf("std(witherr): %e %e\n", av, std);

  free(x);
  return 0;
}
