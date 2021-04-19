#include <stdlib.h>

void rk2(double* xn, int n, double t, double dt, 
         void (*f)(double*, double, double*), double* res)
{
  double *k1 = malloc(n*sizeof(double)); 
  f(xn, t, k1);
  double *xhalf = malloc(n*sizeof(double)); 
  for(int i=0; i<n; ++i) xhalf[i] = xn[i] + k1[i]*dt/2;
  double *k2 = malloc(n*sizeof(double)); 
  f(xhalf, t+dt/2, k2);

  for(int i=0; i<n; ++i) res[i] = xn[i] + dt*k2[i];

  free(k2);
  free(xhalf);
  free(k1);
}

