// to compile: gcc -o executable_name quiz.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printmat(double *m, int N) // print out matrix m which has dimension NxN
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%lf ",m[i*N+j]);
        }
        printf("\n");
    }
}

double trace(double *m, int N) // sum of elements along the diagonal of matrix m (NxN)
{
    double sum = 0;
    for(int i=0; i<N; i++) sum += m[i*N+i];
    return sum;
}

int main(int argc, char *argv[])
{
    double A[4] = {1,2,3,4};
    printmat(A,2);
    printf("Trace: %lf",trace(A,2));
    return 0;
}