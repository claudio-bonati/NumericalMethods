#include<math.h>  // <---- needed for mathematical operations less trivial than basic arithmetical ones
#include<stdio.h>
#include<stdlib.h>

// global variable: its scope is the whole program
// from a softwer engeneering point of view it is better to avoid global 
// variable, since they introduce hidden dependences between functions
// but sometimes are useful
int globint;

void func_on_glob(void)
  {
  globint+=1;
  // equivalent to globint=globint+1;
  }

int product(int a, int b)
  {
  int ris;

  ris=a*b;
  return ris;
  }

// a and b passed by value (not by reference) so this swap does not work!
void swap1(int a, int b)
  {
  int c;
   
  c=a;
  a=b;
  b=c;

  printf("inside swap1: a=%d\n",a);
  printf("inside swap1: b=%d\n",b);
  }


// a and b passed by reference (i.e. with a pointer to their memory location)
void swap2(int *a, int *b)
  {
  int c;
   
  c=*a;
  *a=*b;
  *b=c;

  printf("inside swap2: a=%d\n",*a);
  printf("inside swap2: b=%d\n",*b);
  }


// recursive function
int factorial(int n)
  {
  if(n==1)
    {
    return 1;
    }
  else
    {
    return n*factorial(n-1);
    }
  }


// main
int main(void)
    {
    int a, b, c;
    double x;

    a=2;
    b=3;

    c=product(a,b);

    printf("function 1: product\n");
    printf("%d*%d=%d\n", a,b,c);

    printf("library function: sin\n"); 
    x=2.0;
    printf("%lf\n", sin(x)); // function declared in math.h
                             // in compilation requires -lm to
                             // link to the mathematical library

    printf("\n");

    printf("function 2: wrong swap\n");
    printf("initial a=%d\n",a);
    printf("initial b=%d\n",b);
    swap1(a,b);
    printf("final a=%d\n",a);
    printf("final b=%d\n",b);

    printf("\n");

    printf("function 3: swap\n");
    printf("initial a=%d\n",a);
    printf("initial b=%d\n",b);
    swap2(&a,&b);   
    printf("final a=%d\n",a);
    printf("final b=%d\n",b);

    printf("\n");

    // pay attention that too many recursions could result in a stack overflow 
    // e.g. for a cluster update in the low temperature on large lattices
    printf("function 4: recursive factorial\n");
    c=factorial(5);
    printf("5!=%d\n",c);
    printf("\n");

    // use of global variables
    globint=1;
    printf("globint=%d\n", globint);
    func_on_glob();
    printf("globint=%d\n", globint);

    return EXIT_SUCCESS;
    }


