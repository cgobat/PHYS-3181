#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double gaussrnd()
{
  double u1 = drand48();
  double u2 = drand48();
  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

int main()
{
    for(int i=0;i<100000;i++) printf("%lf\n",gaussrnd());
    return 0;
}