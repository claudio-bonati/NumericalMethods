#include<math.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>

//-------------------------------
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
    {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = (uint32_t) ( ((oldstate >> 18u) ^ oldstate) >> 27u );
    uint32_t rot = (uint32_t) ( oldstate >> 59u );
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
    {
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
    }
//-----------------------------------


//----------------my wrapper for pcg32

// random number internal state
pcg32_random_t pcg32_random_state;

// initialization
void myrand_init(unsigned long int initstate, unsigned long int initseq)
  {
  pcg32_srandom_r(&pcg32_random_state, (uint64_t) initstate, (uint64_t) initseq);
  }

// number in [0,1)
double myrand(void)
  {
  return (double) pcg32_random_r(&pcg32_random_state)/(pow(2.0, 32.0));
  }
//-----------------


int main(int argc, char **argv)
    {
    unsigned long int seed1, seed2; 
    int iterations, i;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s seed1 seed2 iterations\n\n", argv[0]);
      fprintf(stdout, "  seed1, seed2 = seeds for the random numer generator\n");
      fprintf(stdout, "  iterations = number of random numers to be produced\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      seed1=(unsigned long)atol(argv[1]);  // explicit cast to unsigned long to avoid warnings
      seed2=(unsigned long)atol(argv[2]);
      iterations=atoi(argv[3]);
      }

    myrand_init(seed1, seed2);

    for(i=0; i<iterations; i++)
       {
       printf("%lf\n", myrand());
       }

    return EXIT_SUCCESS;
    }


