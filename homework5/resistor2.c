#include <stdio.h>
#include <stdlib.h>
#include "linalgebra.h"

struct resistor
{
  int v1, v2;
  double R;
};

double compute_resistance(int nr, int nv, struct resistor* lr, 
                          int in, int out)
{
  int neq = nr+nv-1;
  double *m = malloc(neq*neq*sizeof(double));
  double *b = malloc(neq*sizeof(double));
  double *x = malloc(neq*sizeof(double));

  for(int i=0; i<neq*neq; i++) m[i] = 0.0;

  // i=0..nr-1  m[i,i] = -Ri
  // i=0..nr-1  m[i,j] = +1 j=nr+inlead_i
  // i=0..nr-1  m[i,j] = -1 j=nr+outlead_i
  for(int i=0; i<nr; ++i)
  {
    m[i*neq+i] = -lr[i].R;
    int j=nr + lr[i].v1;
    if( j != nr+nv-1 ) m[i*neq + j] = +1;
    j = nr + lr[i].v2;
    if( j != nr+nv-1 ) m[i*neq + j] = -1;
  }

  // i=0..nr-1  j = nr + inlead_i m[j,i]=+1
  // i=0..nr-1  j = nr + outlead_i m[j,i]=-1
  // j != nr + nv -1
  for(int i=0; i<nr; ++i)
  {
    int j = nr + lr[i].v1;
    if( j != nr+nv-1 ) m[j*neq+i] = +1;
    j = nr + lr[i].v2;
    if( j != nr+nv-1 ) m[j*neq+i] = -1;
  }

  //printmat(m, neq);
    
  for(int i=0; i<neq; ++i) b[i] = 0;
  int i = nr + in;
  if(i != nr + nv -1) b[i] = +1;
  i = nr + out;
  if(i != nr + nv -1) b[i] = -1;

  gausselim(m, neq, b);
  bks(m, neq, b, x);

  double vin, vout;
  if (in  == nv-1) vin  = 0; else vin  = x[nr+in];
  if (out == nv-1) vout = 0; else vout = x[nr+out];
  double res = vin - vout;

  free(m);
  free(b);
  free(x);

  return res;
}

int main(int argc, char** argv)
{
  int N = atoi(argv[1]);
  int nr = 3*N-2;
  int nv = 2*N;

  printf("N: %d nr: %d nv: %d\n", N, nr, nv);

  struct resistor *lr = malloc(nr*sizeof(struct resistor));

  // bottom resistors
  for(int i=0; i<N-1; ++i) 
  {
    lr[i].v1 = i;
    lr[i].v2 = i+1;
    lr[i].R = 1.0;
  }

  // middle resistors
  for(int i=N-1; i<2*N-1; ++i)
  {
    lr[i].v1 = i-(N-1);
    lr[i].v2 = i-(N-1)+N;
    lr[i].R = 1.0;
  }

  // top resistors
  for(int i=2*N-1; i<3*N-2; ++i)
  {
    lr[i].v1 = i - (N-1);
    lr[i].v2 = i - (N-1)+1;
    lr[i].R = 1.0;
  }

  //for(int i=0; i<nr; ++i) printf("R[%d]: %f  %d %d\n", i, lr[i].R,
  // lr[i].v1, lr[i].v2);

  double r = compute_resistance(nr, nv, lr, 0, N);
  printf("The equivalent resistance is: %e\n", r);

  free(lr);

  return 0;
}
