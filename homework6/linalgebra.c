#include <stdio.h>
#include <math.h>

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

void swap(double* a, double* b)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

void swaprow(double* a, int N, int r1, int r2)
{
  for(int c=0; c<N; ++c) swap(&a[r1*N+c], &a[r2*N+c]);
}

void gausselim(double* a, int N, double* b)
{
  for(int c=0; c<N-1; ++c)
  {
    // scan for the largest pivot
    double piv = fabs(a[c*N+c]); int rm = c;
    for(int r=c+1; r<N; ++r) if( piv < fabs(a[r*N+c]) )
    {
      piv = fabs(a[r*N+c]);
      rm = r;
    }
    if(rm != c) 
    {
      swaprow(a, N, c, rm);
      swap(&b[c], &b[rm]);
    }

    for(int r=c+1; r<N; ++r)
    {
      double fact = a[r*N+c]/a[c*N+c];
      for(int k=c; k<N; ++k) a[r*N+k] -= fact*a[c*N+k];
      b[r] -= fact*b[c];
    }
  }
}

void bks(double* a, int N, double* b, double* x)
{
  for(int i=N-1; i>-1; i--)
  {
    double tmp = b[i];
    for(int j=i+1; j<N; ++j) tmp -= a[i*N+j]*x[j];
    x[i] = tmp/a[i*N+i];
  }
}

void solve(double* a, int N, double* b, double* x)
{
  gausselim(a, N, b);
  bks(a, N, b, x);
}
