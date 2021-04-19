#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double p(double x)
{
  //return exp(-x*x)+2*exp(-(x-5)*(x-5)/.5);

  return exp(-x*x);
}

double delta = 0.5;

void update(double *x)
{
  double step = 2*drand48()-1;
  double nx = (*x) + delta*step;

  double r = p(nx)/p(*x);

  if(r>1) *x=nx;
  if(r<1)
  {
    double rnd = drand48();
    if(rnd < r) *x=nx;
  }
}

int main(int argc, char** argv)
{
  int N = atoi(argv[1]);
  double x = 10;

  for(int i=0; i<N; ++i)
  {
    update(&x);
    printf("%e\n", x);
  }

  return 0;
}
