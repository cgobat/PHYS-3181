#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void gauss_seidel(double* phi, int nphi, int* ne, double exitcrit, int maxiter)
{
  for(int i=0; i<maxiter; ++i)
  {
    double err = 0;
    double nrm = 0;
    for(int j=0; j<nphi; ++j)
    {
      double tmp = 0;
      for(int k=0; k<4; ++k) 
      {
//         printf("k,j: %d %d -> ne[]: %d\n", k, j, ne[4*j+k]); fflush(stdout);
         tmp += phi[ne[4*j+k]];
      }
      err += pow(phi[j]-tmp/4,2);
      nrm += pow(tmp/4, 2);
      phi[j] = tmp/4;
    }
    err = sqrt(err/nrm);
    if(err < exitcrit) break;
  } 
}

struct point {
  int idx;
  int ne[4];
  double pos[2];
};

void setup(int N, int *nphi, int *nb, int** ne, double** pos)
{
  double dx=3.0/(N-1);

  struct point *pts = malloc(N*N*sizeof(struct point));

  for(int i=0; i<N; ++i)
  for(int j=0; j<N; ++j)
  {
    int k = i*N+j;
    pts[k].idx = k;
    pts[k].ne[0] = (i+1)*N+j;
    pts[k].ne[1] = (i-1)*N+j;
    pts[k].ne[2] = i*N+j+1;
    pts[k].ne[3] = i*N+j-1;

    pts[k].pos[0] = i*dx-1.5;
    pts[k].pos[1] = j*dx-1.5;
  }

  double eps = 1e-10;
  for(int i=0; i<N*N; ++i)
  {
    if( (fabs(pts[i].pos[0]) > 1.5-eps) || (fabs(pts[i].pos[1]) > 1.5-eps) )
      pts[i].idx = -1;
    if( (fabs(pts[i].pos[0]) < 0.5+eps) && (fabs(pts[i].pos[1]) < 0.5+eps) )
      pts[i].idx = -2;
  }

  *nphi = 0;
  for(int i=0; i<N*N; ++i) if(pts[i].idx >= 0) 
  {
    pts[i].idx = *nphi;
    (*nphi)++;
  }

  *ne = malloc(4*(*nphi)*sizeof(int));
  *pos = malloc(2*(*nphi)*sizeof(double));

  for(int i=0; i<N*N; ++i) if(pts[i].idx >= 0)
  {
    int offset = pts[i].idx;
    (*ne)[4*offset+0] = pts[pts[i].ne[0]].idx; 
    (*ne)[4*offset+1] = pts[pts[i].ne[1]].idx; 
    (*ne)[4*offset+2] = pts[pts[i].ne[2]].idx; 
    (*ne)[4*offset+3] = pts[pts[i].ne[3]].idx;

    (*pos)[2*offset+0] = pts[i].pos[0]; 
    (*pos)[2*offset+1] = pts[i].pos[1]; 
  } 
  
  *nb = 2;

  free(pts);
}

void print_gradient(int nphi, double* phi, int* ne, double* pos, double a)
{
  for(int i=0; i<nphi; ++i)
  {
    double gx = 0.5*(phi[ne[4*i+0]]-phi[ne[4*i+1]])/a;
    double gy = 0.5*(phi[ne[4*i+2]]-phi[ne[4*i+3]])/a;
    printf("GRAD: %e %e %e %e\n", pos[2*i+0], pos[2*i+1], gx, gy);
  }
}


int main(int argc, char** argv)
{
  int N = atoi(argv[1]);

  int nphi, nb;
  int *ne;
  double *pos;
  setup(N, &nphi, &nb, &ne, &pos);
//printf("nphi: %d nb: %d\n", nphi, nb);
//for(int i=0; i<nphi; ++i) printf("pos[%d]: %e %e\n", i, pos[2*i+0], pos[2*i+1]);

  // phi[i] -> *(phi+i)
  double *data = malloc((nb+nphi)*sizeof(double));
  // phi -> data+nb
  double *phi = data+nb;
  phi[-1] = 100;
  phi[-2] = 0;

  for(int i=0; i<nphi; ++i) phi[i] = 0;
    
  gauss_seidel(phi, nphi, ne, 1e-12, 10000);

  for(int i=0; i<nphi; ++i) printf("PHI: %e %e %e\n", pos[2*i+0], pos[2*i+1], phi[i]);
  
  double a = 3.0/(N-1);
  print_gradient(nphi, phi, ne, pos, a);

  double xmax=0, xmin=1.5;
  int imax,  imin;

  for(int i=0; i<nphi; ++i) if( pos[2*i+0]==pos[2*i+1] )
  {
    double x = pos[2*i+0] ;
    if(x > xmax) { xmax = x; imax = i; }
    if(x>0 && (x<xmin)) { xmin = x; imin = i;}
  }

  double fmax=0, fmin=0;

  for(int i=imin; pos[2*i+0] > -xmin+1e-10 ; i=ne[4*i+1])
  {
    double ey = -(phi[ne[4*i+2]] - phi[ne[4*i+3]])/(2*a);
    fmin += ey*a;
  }

  for(int i=imax; pos[2*i+0] > -xmax+1e-10; i=ne[4*i+1])
  {
    double ey = -(phi[ne[4*i+2]] - phi[ne[4*i+3]])/(2*a);
    fmax += ey*a;
  }

  fmax *= 4;
  fmin *= 4;

  printf("flux(min,max): %.15e %.15e\n", fmin, fmax); 

  free(data);
  free(ne);
  free(pos);

  return 0;
}
