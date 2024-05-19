#include<math.h>
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

    xc='a';
    xi=0;
    xl=0;
    xf=0.0;
    xd=0.0;

    // this is an error
    //xi=0.0;
    //
    // the correct way is using a cast
    //xi=(int) 0.0; 

    xi=xi+one;
    xl=xl+(long)one; // this cast is not really necessary
    xf=(float)one;   // this cast is not really necessary
    xd=1.0;

    // this is an error, since "one" was defined constant
    //one=2;  
    
    printf("%c\n", xc);
    printf("%d\n", xi);
    printf("%ld\n", xl);
    printf("%f\n", xf);
    printf("%lf  %.15lf  %.10e  %.10g\n", xd, xd, xd, xd);

    pxi=&xi;       // pxi point to the location where xi is stored (the "reference" of xi)
                   // *pxi is the content of the memory location pointed by pxi
    *pxi=(*pxi)*2;  
    printf("%d\n",xi);

    return EXIT_SUCCESS;
    }


