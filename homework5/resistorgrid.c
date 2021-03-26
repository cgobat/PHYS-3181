#include <stdio.h>
#include <stdlib.h>
#include "linalgebra.h"

#define RESISTANCE 1.0

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

    for(int i=0; i<neq*neq; i++) m[i] = 0.0;

    for(int i=0; i<nr; ++i)
    {
        m[i*neq+i] = -lr[i].R;
        int j=nr + lr[i].v1;
        if( j != nr+nv-1 ) m[i*neq + j] = +1;
        j = nr + lr[i].v2;
        if( j != nr+nv-1 ) m[i*neq + j] = -1;
    }

    for(int i=0; i<nr; ++i)
    {
        int j = nr + lr[i].v1;
        if( j != nr+nv-1 ) m[j*neq+i] = +1;
        j = nr + lr[i].v2;
        if( j != nr+nv-1 ) m[j*neq+i] = -1;
    }

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

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int x1 = atoi(argv[2]);
    int y1 = atoi(argv[3]);
    int x2 = atoi(argv[4]);
    int y2 = atoi(argv[5]);

    int nr = 2*N*N;
    int nv = N*N;
    int in = x1*N+y1;
    int out = x2*N+y2;
    struct resistor lr[nr];

    for(int i=0; i<nr; i++) lr[i].R = RESISTANCE;

    for(int i = 0; i<N; i++)        // These loops iterate over columns and rows, with the first half of the lr array dedicated to the horiztonal resistors
    {                               // and the second half dedicated to the vertical ones. Although order does not matter as long as all vertices are properly connected,
        for(int j = 0; j<N; j++)    // I have written this enumeration algorithm based on the labels on the diagrams in the assignment. This can be verified by enabling
        {                           // the printout in line 90.
            // horizontals
            lr[i+N*j].v1 = j*N + i;
            lr[i+N*j].v2 = (i==N-1) ? j*N : j*N + i+1; // if at the end of the row, wrap. Otherwise it is (x2,y2)=(x1+1,y1)
            // verticals
            lr[N*N+i+N*j].v1 = j*N + i;
            lr[N*N+i+N*j].v2 = (j==N-1) ? i : (j+1)*N + i; // similar wrapping idea as above, except for the top of each column
        }
    }

    // comment/uncomment the line below to toggle output of the circuit setup itself:
    //for(int i=0; i<nr; i++) printf("Resistor #%d: %.2f Ohm(s), vertices %d and %d\n", i, lr[i].R, lr[i].v1, lr[i].v2);

    double r = compute_resistance(nr, nv, lr, in, out);
    printf("R_eq: %f\n", r);

    return 0;
}