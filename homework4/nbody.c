#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NBODY 3
#define DIM 2

double mass[NBODY];

void compute_force(double *ri, double *rj, double mi, double mj, double *fij)
{
  double G=4*M_PI*M_PI;
  double rij[DIM];
  for(int a=0; a<DIM; ++a) rij[a] = rj[a]-ri[a];
  double rsq=0;
  for(int a=0; a<DIM; ++a) rsq += rij[a]*rij[a];
  double r = sqrt(rsq);
  for(int a=0; a<DIM; ++a) fij[a] = (G*mi*mj/(r*r*r))*rij[a];
}

void f(double *x, double t, double* res)
{
  for(int i=0; i<NBODY; ++i)
  {
    for(int a=0; a<DIM; ++a) res[i*2*DIM+a] = x[i*2*DIM+a+DIM];

    double totalforce[DIM];
    for(int a=0; a<DIM; ++a) totalforce[a] = 0;
   
    for(int j=0; j<NBODY; ++j) if(i!=j)
    {
      double force[DIM];
      compute_force(&x[i*2*DIM], &x[j*2*DIM], mass[i], mass[j], force);
      for(int a=0; a<DIM; ++a) totalforce[a] += force[a];
    }

    for(int a=0; a<DIM; ++a) res[i*2*DIM+a+DIM] = totalforce[a]/mass[i];

  }
}


void rk2(double* xn, double t, double dt, double* res)
{
  double k1[2*NBODY*DIM]; f(xn, t, k1);
  double xhalf[2*NBODY*DIM]; 
  for(int i=0; i<2*NBODY*DIM; ++i) xhalf[i] = xn[i] + k1[i]*dt/2;
  double k2[2*NBODY*DIM]; f(xhalf, t+dt/2, k2);

  for(int i=0; i<2*NBODY*DIM; ++i) res[i] = xn[i] + dt*k2[i];
}

int main(int argc, char** argv)
{
  for(int i = 1; i<argc; i++) printf("%s ",argv[i]);
  double dt = atof(argv[1]);
  double tmax = atof(argv[2]);
  int nstep = tmax/dt;
  double masses[NBODY];
  double positions[DIM*NBODY];
  double velocities[DIM*NBODY];
  for(int i=0;i<3;i++) masses[i] = atof(argv[3+i*5]);
  for(int i = 0; i<DIM*NBODY; i+=DIM)
  {
    for(int j = 0; j<DIM; j++)
    {
      positions[i+j] = atof(argv[4+j+i*5]);
      velocities[i+j] = atof(argv[4+DIM+j+i*5]);
    }
  }
  printf("foo\n");
  for(int i=0; i<NBODY*DIM; ++i) printf("%.3e ", positions[i]);
  // double m1 = atof(argv[3]);
  // double x1 = atof(argv[4]);
  // double y1 = atof(argv[5]);
  // double vx1 = atof(argv[6]);
  // double vy1 = atof(argv[7]);
  // double m2 = atof(argv[8]);
  // double x2 = atof(argv[9]);
  // double y2 = atof(argv[10]);
  // double vx2 = atof(argv[11]);
  // double vy2 = atof(argv[12]);
  // double m3 = atof(argv[13]);
  // double x3 = atof(argv[14]);
  // double y3 = atof(argv[15]);
  // double vx3 = atof(argv[16]);
  // double vy3 = atof(argv[17]);
  // double state[2*NBODY*DIM] = {x1, vx1, y1, vy1, x2, vx2, y2, vy2, x3, vx3, y3, vy3};

printf("%.3e ", 0.0);
for(int i=0; i<2*NBODY*DIM; ++i) printf("%.3e ", state[i]);
printf("\n");

for(int i=0; i<nstep; ++i)
{
  double result[2*NBODY*DIM];
  rk2(state, i*dt, dt, result);
  printf("%.3e ", i*dt+dt);
  for(int j=0; j<2*NBODY*DIM; ++j) printf("%.3e ", result[j]);
  printf("\n");
  for(int j=0; j<2*NBODY*DIM; ++j) state[j] = result[j];
}

  return 0;
}
