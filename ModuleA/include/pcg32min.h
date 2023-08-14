#ifndef PCG32MIN_H
#define PCG32MIN_H

#include<math.h>
#include<stdint.h>

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

extern pcg32_random_t pcg32_random_state;

uint32_t pcg32_random_r(pcg32_random_t* rng);
void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq);


#endif

