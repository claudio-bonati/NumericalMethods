#ifndef RANDOM_H
#define RANDOM_H

// initialize the random number generator
void myrand_init(unsigned long int initstate, unsigned long int initseq);

// return a random number in [0,1)
double myrand(void);

#endif
