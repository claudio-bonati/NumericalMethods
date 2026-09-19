#include<complex.h> 
#include<stdio.h>
#include<stdlib.h>

// this code needs -std=c99 to be added to the compilation string
// to properly compile, although it is likely that some compiler
// can understand this also without explicitly specifying it 

// main
int main(void)
    {
    double complex x, y,z;

    x=2.0+0.0*I;
    y=0.0+1.0*I;
    z=0.0+1.0*I;

    printf("%lf %lf\n", creal(x*y*z), cimag(x*y*z));

    return EXIT_SUCCESS;
    }


