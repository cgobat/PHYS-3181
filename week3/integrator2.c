#include <stdlib.h>
#include <stdio.h>

void f(double *x, double t, double* res)
{
  res[0] = x[2];
  res[1] = x[3];
  res[2] = -2.*x[0]*x[1]*x[1];
  res[3] = -2.*x[0]*x[0]*x[1];
}


void rk2(double* xn, double t, double dt, double* res)
{
  double k1[4]; f(xn, t, k1);
  double xhalf[4]; 
  for(int i=0; i<4; ++i) xhalf[i] = xn[i] + k1[i]*dt/2;
  xhalf[1] = xn[1] + k1[1]*dt/2;
  double k2[4]; f(xhalf, t+dt/2, k2);

  for(int i=0; i<4; ++i) res[i] = xn[i] + dt*k2[i];
}

int main(int argc, char** argv)
{
  double dt = atof(argv[1]);
  int nstep = atoi(argv[2]);
//printf("argc: %d dt: %f nstep: %d\n", argc, dt, nstep);

  double x[4] = {0,0,1,2};

  printf("%.5lf %.5lf %.5lf %.5lf %.5lf\n", 0.0, x[0], x[1], x[2], x[3]);

  for(int i=0; i<nstep; ++i)
  {
    double xf[4];
    rk2(x, i*dt, dt, xf);
    //double xnp1 = euler(x, i*dt, dt);
    printf("%.5lf %.5lf %.5lf %.5lf %.5lf\n", i*dt+dt, xf[0], xf[1], xf[2], xf[3]);
    for(int j=0; j<4; ++j) x[j] = xf[j];
  }

  return 0;
}
