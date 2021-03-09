#include <stdio.h>
#include <stdlib.h>
#include "linalgebra.h"

struct resistor
{
    int v1, v2;
    double R;
};

double compute_resistance(int nr, int nv, struct resistor* lr, int in, int out)
{
	int neq = nr+nv-1;
	double *m = malloc(neq*neq*sizeof(double));
	double *b = malloc(neq*sizeof(double));
	double *x = malloc(neq*sizeof(double));

	for(int i=0; i<neq*neq; i++) m[i] = 0;

	for(int i=0; i<nr; i++)
	{
		m[i*neq+i] = -lr[i].R;
		if((lr[i].v2) != (nv-1)) m[i*neq+lr[i].v2+nr] = 1.0;
		if((lr[i].v1) != (nv-1)) m[i*neq+lr[i].v1+nr] = -1.0;
	}

	for(int i=0; i<nr; i++)
	{
		if(lr[i].v1 != (nv-1)) m[(nr+lr[i].v1)*neq + i] = -1.0;
		if(lr[i].v2 != (nv-1)) m[(nr+lr[i].v2)*neq + i] = 1.0;
	}

	for(int i=0; i<neq; ++i) b[i] = 0;
	if(in != (nv-1)) b[nr+in] = -1.0;
	if(out != (nv-1)) b[nr+out] = +1.0;

	gausselim(m,neq,b);
	bks(m,neq,b,x);

	//printmat(m,neq);

	double vin = (in == nv-1)? 0 : x[nr+in];
	double vout = (out == nv-1)? 0 : x[nr+out];
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

	struct resistor lr[] = {0,1,1.0,  1,2,2.0,  3,2,3.0,  0,3,4.0,  1,3,5.0};

	for(int i=0; i<nr; ++i)
		printf("R[%d]: %f %d %d\n", i, lr[i].R, lr[i].v1, lr[i].v2);

	double r = compute_resistance(nr, nv, lr, 0, 2);
	printf("Equivalent resistance is: %e\n",r);

	return 0;
}