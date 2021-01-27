// to compile: gcc quiz1.c -o quiz1

#include <stdio.h>

int main()
{
    double input;
    printf("Enter the number to square: ");
    scanf("%lf", &input);
    
    printf("%lf squared is %lf",input,input*input);

    return 0;
}