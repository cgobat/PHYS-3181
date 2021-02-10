#include <stdio.h>

int sumton(int n)
{
    int sum = 0;
    for(int i = 0; i<n+1; i++) sum += i;
    return sum;
}

int main()
{
    for(int i=1; i<101; i++)
    {
        printf("%d\t%d\n",i,sumton(i));
    }
    return 0;
}