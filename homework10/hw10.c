#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 30. // dimension of main square
#define RADIUS 10. // radius of center hole

struct point // using a data structure to combine all of each point's info will make setup much easier
{
    int idx; // index in the array
    double x,y; // physical x and y positions
    int ne[4]; // indices of neighboring points
};

void initialize(int N, int *nphi, int *ne, double *pos)
{
    struct point* p = malloc(N*N*sizeof(struct point));

    double a = SIZE /(N-1); // spacing

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
        int idx = i*N+j;
        p[idx].idx = i*N+j;
        p[idx].x = j*a - SIZE/2.0 ;
        p[idx].y = i*a - SIZE/2.0 ;
        p[idx].ne[0] = i*N+j+1;
        p[idx].ne[1] = i*N+j-1;
        p[idx].ne[2] = (i+1)*N+j;
        p[idx].ne[3] = (i-1)*N+j;
    }

    double eps = 1e-10;
    for(int i=0; i<N*N; i++)
    {
        if(fabs(p[i].x) > SIZE/2. - eps) p[i].idx = -1;
        if(fabs(p[i].y) > SIZE/2. - eps) p[i].idx = -1;
        if(pow(fabs(p[i].x),2) + pow(fabs(p[i].y),2) < RADIUS*RADIUS) p[i].idx=-2;
    }

    int newidx=0;
    for(int i=0; i<N*N; i++) if(p[i].idx >= 0) p[i].idx=newidx++;
    *nphi = newidx;

    for(int i=0; i<N*N; i++) if(p[i].idx >= 0)
    {
        int idx = p[i].idx;
        pos[2*idx+0] = p[i].x;
        pos[2*idx+1] = p[i].y;
        ne[4*idx+0] = p[p[i].ne[0]].idx;
        ne[4*idx+1] = p[p[i].ne[1]].idx;
        ne[4*idx+2] = p[p[i].ne[2]].idx;
        ne[4*idx+3] = p[p[i].ne[3]].idx;
    }

    free(p);
}

int gauss_seidel(int nphi, double* phi, int *ne, double tolerance)
{
    int itercount = 0;
    double err=0, nrm=0;
    do {    
        for(int i=0; i<nphi; i++)
        { 
            double tmp = (phi[ne[4*i+0]]+ phi[ne[4*i+1]] + phi[ne[4*i+2]] + phi[ne[4*i+3]])/4;
            err += (tmp-phi[i])*(tmp-phi[i]);
            nrm += tmp*tmp;
            phi[i] = tmp;
        }
	    itercount++;
    } while (sqrt(err/nrm) > tolerance);
    // printf("Error: %e,\tnrm: %e\n", sqrt(err), sqrt(nrm));
    return itercount;
}

int main(int argc, char *argv[])
{
    int iter;
    int N = atoi(argv[1]); // resolution (points per side of square)
    double delta = atof(argv[2]); // tolerance (desired precision)

    int* ne = malloc(4*N*N*sizeof(int)); // neighbors
    double* pos = malloc(2*N*N*sizeof(double)); // positions

    int nphi;
    initialize(N, &nphi, ne, pos);

    // for(int i=0; i<nphi; i++) // verify neighbors are correct
    //     printf("ne[%d]: %d %d %d %d\n", i, ne[4*i+0], ne[4*i+1], ne[4*i+2], ne[4*i+3]);

    double* phi = malloc((nphi+2)*sizeof(double));
    phi[0] = 500; // boundary condition, metal pipe
    phi[1] = 300; // boundary condition, room temp
    for(int i=0; i<nphi; i++) phi[2+i] = 500;

    iter = gauss_seidel(nphi, phi+2, ne, delta); // note: phi+2 is like &(phi[2]) -- same effect as passing phi[2:] in python
    
    for(int i=0; i<nphi; i++)
        printf("Result @ (%e, %e): %e\n", pos[2*i+0], pos[2*i+1],phi[2+i]);
    
    printf("Desired precision attained in %d iterations.\n",iter);

    double a = SIZE/(N-1);

    for(int i=0; i<nphi; i++)
    { // compute gradients
        double ex = -(phi[2+ne[4*i+0]] - phi[2+ne[4*i+1]])/(2*a);
        double ey = -(phi[2+ne[4*i+2]] - phi[2+ne[4*i+3]])/(2*a);
    }

    int max_i, min_i;
    double max_x = 0., min_x = 30.;

    for(int i=0; i<nphi; i++) if(pos[2*i+0] == pos[2*i+1]) // x=y
    {
	    double x = pos [2*i+0];
	    if(x>max_x) { max_x=x; max_i = i; }
	    if(x>0 && x<min_x) { min_x=x; min_i = i; }
    }

    double max_f = 0., min_f = 0.;

    for(int i=max_i; pos[2*i+0] > -max_x+a/2; i = ne[4*i+1])
    {
	    double ey = -(phi[2+ne[4*i+2]]-phi[2+ne[4*i+3]]) /(2*a);
	    max_f += a*ey;
    }
    
    for(int i=min_i; pos[2*i+0] > -min_x+a/2; i = ne[4*i+1])
    {
	    double ey = -(phi[2+ne[4*i+2]]-phi[2+ne[4*i+3]]) /(2*a);
	    min_f += a*ey;
    }

    min_f *= 4;
    max_f *= 4;
    
    printf("Min. flux: %e\tMax. flux: %e\n", min_f , max_f);

    free(ne);
    free(pos);
    free(phi); 
    return 0;
}
