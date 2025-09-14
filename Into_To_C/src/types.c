#include<stdio.h>
#include<stdlib.h>

// main
int main(void)
    {
    char xc;
    int xi;      // integer
    long int xl; // long integer
    float xf;    // single precision floating point
    double xd;   // double precision floating point
    const int one=1; // const integer
    int *pxi;    // pointer to integer

    // Since the C99 standad complex variables are also available, e.g.
    // 
    // double complex x=0.9+0.3*I;
    //
    // To use complex numbers we need to include the <complex.h> header and
    // pass to the compiler the flag -std=c99, to activate the functionalities
    // of the C99 standard
    //
    // see types_complex.c for an explicit example

    xc='a';
    xi=0;
    xl=0;
    xf=0.0;
    xd=0.0;

    // this would be an error (or a warning, triggered by -Wconversion)
    // xi=xd;
    //
    // the correct way is to use a cast
    // xi=(int) xd; 

    xi=xi+one;
    xl=xl+(long)one; // this cast is not really necessary
    xf=(float)one;   // this cast is not really necessary
    xd=100.0;

    xi=2*xi; // this is to be read as: "store in xi the number 2*xi" 
             // NOT as an equation

    // this is an error, since "one" was defined constant
    //one=2;  
    
    printf("1) %c\n", xc);
    printf("2) %d\n", xi);
    printf("3) %ld\n", xl);
    printf("4) %f\n", xf);
    printf("5) %lf  %.15lf  %.10e  %.10g\n", xd, xd, xd, xd);

    pxi=&xi;       // pxi point to the location where xi is stored (the "reference" of xi)
                   // *pxi is the content of the memory location pointed by pxi
                   
    *pxi=(*pxi)*2;  
    printf("6) %d\n",xi);

    return EXIT_SUCCESS;
    }


