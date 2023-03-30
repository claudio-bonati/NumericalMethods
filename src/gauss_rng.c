#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/random.h"

// main
int main (int argc, char **argv)
    {
    long int i, sample;
    double tmp1, tmp2, x, x2, x4, sigma;
    const unsigned long int seed1=time(NULL);
    const unsigned long int seed2=seed1+127;


    if(argc != 2)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s sample\n", argv[0]);
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  test for <x> estimated by using 'sample' Box-Muller draws\n");
      fprintf(stdout, "  test for <x^2> estimated by using 'sample' Box-Muller draws\n");
      fprintf(stdout, "  test for <x^4> estimated by using 'sample' Box-Muller draws\n");

      return EXIT_SUCCESS;
      }
    else
      {
      sample=atol(argv[1]);
      }

    if(sample<0)
      {
      fprintf(stderr, "'sample' mast be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    for(i=0; i<sample/2; i++)
       {
       gauss2(&tmp1, &tmp2);
   
       x+=(tmp1+tmp2);
       x2+=(pow(tmp1,2)+pow(tmp2,2));
       x4+=(pow(tmp1,4)+pow(tmp2,4));
       }
    if(sample % 2 == 1)
      {
      tmp1+=gauss1();

      x+=tmp1;
      x2+=pow(tmp1,2);
      x4+=pow(tmp1,4);
      } 

    x/=(double) sample;
    x2/=(double) sample;
    x4/=(double) sample;
 
    sigma=1/sqrt(sample);
    printf("test <x>: %lf %lf %lf\n", x, sigma, x/sigma);

    sigma=sqrt(3-1)/sqrt(sample);
    printf("test <x^2>: %lf %lf %lf\n", x2-1, sigma, (x2-1)/sigma);

    sigma=sqrt(105-9)/sqrt(sample);
    printf("test <x^4>: %lf %lf %lf\n", x4-3, sigma, (x4-3)/sigma);

    return EXIT_SUCCESS;
    }

