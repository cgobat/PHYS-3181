// Caden Gobat, PHYS 3181
// Code Assignment 4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NBODY 3
#define DIM 2

void compute_force(double *ri, double *rj, double mi, double mj, double *fij)
{
  double G=4*M_PI*M_PI;
  double rij[DIM]; // direction and distance between the two bodies
  for(int a=0; a<DIM; ++a) rij[a] = rj[a]-ri[a];
  double rsq=0; 
  for(int a=0; a<DIM; ++a) rsq += rij[a]*rij[a];
  double r = sqrt(rsq); // AKA the magnitude of the rij vector
  for(int a=0; a<DIM; ++a) fij[a] = (G*mi*mj/(r*r*r))*rij[a];
}

void f(double *x, double t, double* res, double* mass)
{
  for(int i=0; i<NBODY; ++i) // for each body in the system
  {
    for(int a=0; a<DIM; ++a) res[i*2*DIM+a] = x[i*2*DIM+a+DIM];

    double totalforce[DIM]; // net force that body i is subject to
    for(int a=0; a<DIM; ++a) totalforce[a] = 0;
   
    for(int j=0; j<NBODY; ++j) if(i!=j) // for each other body it interacts with
    {
      double force[DIM];
      compute_force(&x[i*2*DIM], &x[j*2*DIM], mass[i], mass[j], force); // compute the force between them
      for(int a=0; a<DIM; ++a) totalforce[a] += force[a]; // add it to the sum of forces that body i feels
    }

    for(int a=0; a<DIM; ++a) res[i*2*DIM+a+DIM] = totalforce[a]/mass[i]; // return the acceleration since this is a derivatives function

  }
}

void rk2(double* xn, double t, double dt, double* res, double* mass)
{
  double k1[2*NBODY*DIM]; f(xn, t, k1, mass);
  double xhalf[2*NBODY*DIM]; 
  for(int i=0; i<2*NBODY*DIM; ++i) xhalf[i] = xn[i] + k1[i]*dt/2;
  double k2[2*NBODY*DIM]; f(xhalf, t+dt/2, k2, mass);

  for(int i=0; i<2*NBODY*DIM; ++i) res[i] = xn[i] + dt*k2[i];
}

int main(int argc, char *argv[])
{
    // the following code block is currently commented out but can be used for debugging the arguments passed:
    // printf("%d/%d args passed: ",argc-1,2+NBODY+2*DIM*NBODY);
    // for(int i = 1; i<argc-1; i++) printf("%s ",argv[i]);
    // printf("\n");


    // begin argument parsing (ordering the arguments as requested in the assignment makes things a little complicated)
    double dt = atof(argv[1]);
    double tmax = atof(argv[2]);
    
    double masses[NBODY];
    double positions[DIM*NBODY];
    double velocities[DIM*NBODY];

    for(int i=0;i<NBODY;i++) masses[i] = atof(argv[3+i*(DIM+NBODY)]); // make an array of just the masses
    
    for(int i = 0; i<DIM*NBODY; i+=DIM)
    {
        for(int j = 0; j<DIM; j++)
        {
            positions[i+j] = atof(argv[4+j+i*(DIM+NBODY)]); // array of positions: for 2D it is [x1, y1, x2, y2, x3, y3, ..., etc.]
            velocities[i+j] = atof(argv[4+DIM+j+i*(DIM+NBODY)]); // array of positions: for 2D it is [vx1, vy1, vx2, vy2, vx3, vy3, ..., etc.]
        }
    }
    
    double state[2*DIM*NBODY]; // now combine positions and velocities into a single state vector
    for(int i=0; i<DIM*NBODY; i++) // this will create (for 2D, N=3) the vector [x1 y1 x2 y2 x3 y3 vx1 vy1 vx2 vy2 vx3 vy3]
    {
        state[i] = positions[i];
        state[i+(DIM*NBODY)] = velocities[i];
    }
    
    // end argument parsing, begin computation and printout

    printf("%.3e ", 0.0);
    for(int i=0; i<2*NBODY*DIM; ++i) printf("%.3e ", state[i]); // print initial conditions (time 0)
    printf("\n");

    for(double t=0; t<tmax; t+=dt)
    {
        double result[2*NBODY*DIM];
        rk2(state, t, dt, result, masses);
        printf("%.3e ", t+dt);
        for(int j=0; j<2*NBODY*DIM; j++) printf("%.3e ", result[j]);
        printf("\n");
        for(int j=0; j<2*NBODY*DIM; j++) state[j] = result[j];
    }
    return 0;
}