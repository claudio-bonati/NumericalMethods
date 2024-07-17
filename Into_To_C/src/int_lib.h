#ifndef INT_LIB_H
#define INT_LIB_H

// compute the integral \int_{xmin}^{xmax) f(x) dx 
// with relative accuracy relacc_target
// using the rectangle method
double int_rel_rec(double xmin, 
                   double xmax, 
                   double (*f)(double), // function pointer
                   double relacc_target,
                   int *err);
 
// compute the integral \int_{xmin}^{xmax) f(x) dx 
// with  accuracy absacc_target
// using the rectangle method
double int_abs_rec(double xmin, 
                   double xmax, 
                   double (*f)(double), // function pointer
                   double absacc_target,
                   int *err);

#endif 
