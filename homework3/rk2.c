#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define GM 4.*M_PI*M_PI // when using units of AU and years

void f(double* x, double* res) // gravitational derivatives function
{
    res[0] = x[2]; // x becomes vx
    res[1] = x[3]; // y becomes vy
    res[2] = GM/(x[0]*x[0] + x[1]*x[1])*-x[0]/sqrt(x[0]*x[0]+x[1]*x[1]); // vx becomes ax
    res[3] = GM/(x[0]*x[0] + x[1]*x[1])*-x[1]/sqrt(x[0]*x[0]+x[1]*x[1]); // vy becomes ay
}

void rk2(double* init, double dt, double* res) // Runge-Kutta integrator function
{
    double k1[4];
    f(init, k1); // perform initial calculation and write to the 4-vector k1

    double half[4]; 
    for(int i=0; i<4; i++)
        half[i] = init[i] + k1[i]*dt/2; // Runge-Kutta half step
    
    double k2[4];
    f(half, k2); // perform second-order calculation using half step and write to the 4-vector k2

    for(int i=0; i<4; ++i)
        res[i] = init[i] + dt*k2[i]; // perform the actual stepping operation
}

int main(int argc, char* argv[])
{
    double x0 = atof(argv[1]); // initial x position
    double y0 = atof(argv[2]); // initial y position
    double vx0 = atof(argv[3]); // initial x velocity
    double vy0 = atof(argv[4]); // initial y velocity
    double dt = atof(argv[5]); // time step size
    double T = atof(argv[6]); // total time to run

    double dts[7] = {0.0001,0.0003,0.001,0.003,0.01,0.03,0.1};
    
    double result[4];

    // for(int i=0; i<7; i++)
    // {
    //     double dt = dts[i];

         double state[4] = {x0,y0,vx0,vy0};

        //printf("%e %e %e %e %e\n", 0., state[0], state[1], state[2], state[3]); // conditions before taking any steps
        
        for (double t = 0; t<T; t+=dt) // only go up to t<T because the output is the state AFTER taking each step
        {
            rk2(state, dt, result); // call the RK step function on the state vector and write output to result vector
            //printf("%lf %lf\n",t/5.,floor(t/5.));
            //if (t/5.-floor(t/5.)<0.01)
              //  printf("%e %e %e %e %e\n", t+dt, result[0], result[1], result[2], result[3]);
            for(int j=0; j<4; ++j) state[j] = result[j];
            if (t/5.-floor(t/5.)<0.01)
                printf("%e %e %e\n", t, 0.5*(state[2]*state[2]+state[3]*state[3]), -GM/sqrt(state[0]*state[0]+state[1]*state[1]));
        }
        //printf("%e %e %e %e\n",state[0],state[1],state[2],state[3]);
        //
    //}
    
    return 0;
}