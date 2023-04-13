#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/random.h"

// main
int main(int argc, char **argv)
    {
    long int i, sample;
    double mu, tmp1, tmp2, aux, x, x2, Z, Z2, sigmax, sigmaZ;
    double ris, err;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 3)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s sample mu\n\n", argv[0]);
      fprintf(stdout, "  sample = number of draws\n");
      fprintf(stdout, "  mu = average of the normal distribution N(mu,1) of which we want to estimate <x>_mu\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  value of <x>_mu obtained by sampling a normal distribution N(0,1)\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input

      sample=atol(argv[1]);
      mu=atof(argv[2]);
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // initialize average values
    x=0.0;
    x2=0.0;
    Z=0.0;
    Z2=0.0;

    // average values are computed by using
    // \int x e^{-(x-mu)^2/2} = \int x e^{-(x-mu)^2/2}e^{x^2/2}e^{-x^2/2} 
    // = (\int e^{-x^2/2}) * <x e^{-(x-mu)^2/2}e^{x^2/2}>_{N(0,1)}
    //
    // \int e^{-(x-mu)^2/2} = \int e^{-(x-mu)^2/2}e^{x^2/2}e^{-x^2/2} 
    // = (\int e^{-x^2/2}) * <e^{-(x-mu)^2/2}e^{x^2/2}>_{N(0,1)}
    //
    // hence <x>_{N(mu,1)}=
    // <x e^{-(x-mu)^2/2}e^{x^2/2}>_{N(0,1)} / <e^{-(x-mu)^2/2}e^{x^2/2}>_{N(0,1)}

    // numerator 
    for(i=0; i<sample/2; i++)
       {
       gauss2(&tmp1, &tmp2);
   
       aux=tmp1*exp(-(tmp1-mu)*(tmp1-mu)*0.5+tmp1*tmp1*0.5);
       x+=aux;
       x2+=aux*aux;

       aux=tmp2*exp(-(tmp2-mu)*(tmp2-mu)*0.5+tmp2*tmp2*0.5);
       x+=aux;
       x2+=aux*aux;
       }
    if(sample % 2 == 1)
      {
      tmp1+=gauss1();

      aux=tmp1*exp(-(tmp1-mu)*(tmp1-mu)*0.5+tmp1*tmp1*0.5);
      x+=aux;
      x2+=aux*aux;
      } 

    x/=(double) sample;
    x2/=(double) sample;
    sigmax=sqrt(x2-x*x)/sqrt((double) sample);

    // denominator
    for(i=0; i<sample/2; i++)
       {
       gauss2(&tmp1, &tmp2);
   
       aux=exp(-(tmp1-mu)*(tmp1-mu)*0.5+tmp1*tmp1*0.5);
       Z+=aux;
       Z2+=aux*aux;

       aux=exp(-(tmp2-mu)*(tmp2-mu)*0.5+tmp2*tmp2*0.5);
       Z+=aux;
       Z2+=aux*aux;
       }
    if(sample % 2 == 1)
      {
      tmp1+=gauss1();

      aux=exp(-(tmp1-mu)*(tmp1-mu)*0.5+tmp1*tmp1*0.5);
      Z+=aux;
      Z2+=aux*aux;
      } 

    Z/=(double) sample;
    Z2/=(double) sample;
    sigmaZ=sqrt(Z2-Z*Z)/sqrt((double) sample);
 
    ris=x/Z;
    err=ris*sqrt( pow(sigmax/x,2.0)+pow(sigmaZ/Z,2.0));
    // we can use error propagation since we used two different random samples

    printf("<x>=%lf %lf\n", ris, err);

    return EXIT_SUCCESS;
    }

