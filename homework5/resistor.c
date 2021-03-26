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

int main()
{
  int nr = 5;
  int nv = 4;

  // lr[0].v1=0; lr[0].v2=1; lr[0].R=1.0; ...
  struct resistor lr[5] = {0,1,1.0,  1,2,2.0, 3,2,3.0, 0,3,4.0,
     1,3,5.0};


  for(int i=0; i<nr; ++i) printf("R[%d]: %f  %d %d\n", i, lr[i].R,
   lr[i].v1, lr[i].v2);

  double r = compute_resistance(nr, nv, lr, 0, 2);
  printf("The equivalent resistance is: %e\n", r);

  return 0;
}
