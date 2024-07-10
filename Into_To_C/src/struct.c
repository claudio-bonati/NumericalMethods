#include<stdio.h>
#include<stdlib.h>

// it is convenient to use structures when several variables often appears together 
// or have some reason to stay together

// struct + typedef 
typedef struct Animal {
  int legs;
  double weight;  
  struct Animal *mother; // this will not be used. I added it only to show that it is possible
                         // to define recursive structures, which are useful to implement
                         // trees and graphs
                         // In the specific case the pointer *mother could be used to point to the 
                         // structure which describes the "mother" or the animal described.
  } Animal;

// the disentangled struct+typedef definition is as follows
//struct Animal {
//  int legs;
//  double weight;  
//  };
//
// at this point to define a variable "dog" you should use
// struct Animal dog;
//
// using 
// typedef struct Animal Animal;
// we create the alias "Animal" for "struct Animal"


double av_w_on_leg(Animal bobby)
  {
  double x;
  x=bobby.weight / bobby.legs;
  return x;
  }


// to avoid copying large structures (remember that function arguments are passed by value)
// it is convenient to use
double av_w_on_leg_large_struct(Animal *bobby)
  {
  double x;
  x=bobby->weight / bobby->legs;  // bobby->weight is equivalent to (*bobby).weight
  return x;
  }

// in case the function does not have to change the elements of "bobby" but only change its
// value it is better to use 
//double av_w_on_leg_large_struct(Animal const * const bobby)
//
// note that using const * const both the pointer and pointee are constant see
// https://en.wikipedia.org/wiki/Const_(computer_programming)


// this function change some element of "bobby", so we can not use const * const
void ten_k_fatter(Animal * bobby)
  {
  bobby->weight+=10.0;
  }


// main
int main(void)
    {
    Animal dog, chicken, lame_cat;

    dog.legs=4;
    dog.weight=50.0;

    chicken.legs=2;
    chicken.weight=5.0;

    lame_cat.legs=3;
    lame_cat.weight=10.0;

    printf("dog=%lf  chicken=%lf  lame_cat=%lf\n", av_w_on_leg(dog), 
                                                   av_w_on_leg(chicken), 
                                                   av_w_on_leg(lame_cat));
    printf("dog=%lf  chicken=%lf  lame_cat=%lf\n", av_w_on_leg_large_struct(&dog), 
                              av_w_on_leg_large_struct(&chicken), 
                              av_w_on_leg_large_struct(&lame_cat)); 

    ten_k_fatter(&chicken);
    printf("fat chicken=%lf\n", av_w_on_leg(chicken));

    return EXIT_SUCCESS;
    }


