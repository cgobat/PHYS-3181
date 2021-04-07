#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double gaussrnd()
{
  double u1 = drand48();
  double u2 = drand48();
  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

void log_dist(double* a, double* b, double* res,  int N)
{
  for(int i=0; i<N; ++i) res[i] = log(a[i]*a[i]+b[i]*b[i]);
}

int main(int argc, char** argv)
{
  int N = atoi(argv[1]);
  double* x = malloc(N*sizeof(double));
  double* y = malloc(N*sizeof(double));
  double* log = malloc(N*sizeof(double));

  for(int i=0; i<N; ++i) x[i] = 1 + 2*gaussrnd();

  for(int i=0; i<N; ++i) y[i] = 4 + 1*gaussrnd();

  log_dist(x, y, log, N);

  for(int j=0; j<N; j++) printf("%lf\n", log[j]);
  free(x);
  free(y);
  free(log);

  return 0;
}