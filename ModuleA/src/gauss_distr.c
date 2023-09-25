#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/random.h"


// test to check if the polar form of Box-Muller is more efficient than the basic one
void test_speed(void)
  {
  long int i, sample=500000000;
  double tmp1, tmp2, diff_sec;
  clock_t time1, time2;

  printf("---------------------\n");
  printf("Test for Box-Muller speed:\n");

  time1=clock();
  for(i=0; i<sample; i++)
     {
     gauss2(&tmp1, &tmp2);
     }
  time2=clock();
  diff_sec=((double) (time2-time1))/CLOCKS_PER_SEC;
  printf("polar form: %.3lf s\n", diff_sec);

  time1=clock();
  for(i=0; i<sample; i++)
     {
     gauss2_basic(&tmp1, &tmp2);
     }
  time2=clock();
  diff_sec=((double) (time2-time1))/CLOCKS_PER_SEC;
  printf("basic form: %.3lf s\n", diff_sec);
  printf("---------------------\n");
  }


// main
int main(int argc, char **argv)
    {
    long int i, sample;
    double tmp1, tmp2, x, x2, x4, sigma;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 2)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s sample\n\n", argv[0]);
      fprintf(stdout, "  sample = number of draws\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  test for <x> estimated by using 'sample' Box-Muller draws\n");
      fprintf(stdout, "  test for <x^2> estimated by using 'sample' Box-Muller draws\n");
      fprintf(stdout, "  test for <x^4> estimated by using 'sample' Box-Muller draws\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      sample=atol(argv[1]);
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // to test which Box-Muller implementation is faster
    //test_speed();

    // initialize average values
    x=0.0;
    x2=0.0;
    x4=0.0;

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
 
    sigma=1.0/sqrt((double)sample);
    printf("<x>-exact=%f ; ", x);
    printf("sigma_th=%f ; ", sigma); // theoretical std
    printf("(<x>-exact)/sigma_th=%f\n", x/sigma);

    sigma=sqrt(3.0-1.0)/sqrt((double) sample);
    printf("<x^2>-exact=%f ; ", x2-1);
    printf("sigma_th=%f ; ", sigma); // theoretical std
    printf("(<x^2)-exact)/sigma_th=%f\n", (x2-1)/sigma);

    sigma=sqrt(105.0-9.0)/sqrt((double) sample);
    printf("<x^4>-exact=%f ; ", x4-3);
    printf("sigma_th=%f ; ", sigma); // theoretical std
    printf("(<x^4>-exact)/sigma_th=%f\n", (x4-3)/sigma);

    return EXIT_SUCCESS;
    }

