#include <stdio.h>

void printmat(double* m, int N)
{
  for(int i=0; i<N; ++i)
  {
    for(int j=0; j<N; ++j) printf("%5.3f ", m[i*N+j]);
    printf("\n");
  }
}

double trace(double* m, int N)
{
  double tr = 0.0;
  for(int i=0; i<N; ++i) tr += m[i*N+i];
  return tr;
}

void gausselim(double* a, int N, double* b)
{
  for(int c=0; c<N-1; ++c)
  {
    for(int r=c+1; r<N; ++r)
    {
      double fact = a[r*N+c]/a[c*N+c];
      for(int k=c; k<N; ++k) a[r*N+k] -= fact*a[c*N+k];
      b[r] -= fact*b[c];
    }
  }
}
