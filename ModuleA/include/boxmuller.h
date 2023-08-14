#ifndef BOXMULLER_H
#define BOXMULLER_H

// single normal random number 
// polar form
double gauss1(void);

// couple of independent normal random numbers
// polar form
void gauss2(double * restrict ris1, double * restrict ris2);

// couple of independent normal random numbers
// basic form
void gauss2_basic(double * restrict ris1, double * restrict ris2);

#endif
