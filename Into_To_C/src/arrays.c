#include<stdio.h>
#include<stdlib.h>
#include<string.h> // <---- needed to use string operations

#define LENGTH 10  // macro: at compile time LENGTH is 
                   // replaced in the code by its numerical value
// if LENGTH is too large one gets in execution
// Segmentation fault (core dumped)
// and dynamic memory allocation is required
// this happens before the RAM is full


void times2(int x[LENGTH]) // also "void times2(int *x)" would work
   {
   int i;
 
   for(i=0; i<LENGTH; i++)
      {
      x[i]*=2;
      }
   }

// void times2wrong(int x[LENGTH])
//
// could also be used insted of
//
// void times2wrong(int *x)
//
// and in fact it would be better: the compiler could detect the problem in
// the function and complain
void times2wrong(int *x)
   {
   int i;
 
   for(i=0; i<LENGTH+10; i++) // i reaches values > LENGTH. Problem!
      {
      x[i]*=2;
      }
   }


void times2matrix(int x[LENGTH][LENGTH]) // times2matrix(int **x) can produce warnings
   {
   int i, j;
 
   for(i=0; i<LENGTH; i++)
      {
      for(j=0; j<LENGTH; j++)
         {
         x[i][j]*=2;
         }
      }
   }


// main
int main(void)
    {
    int i, j;
    int v[LENGTH];
    int matrix[LENGTH][LENGTH];
    char name[20];

    // -------------------- vector 

    for(i=0; i<LENGTH; i++)
       {
       v[i]=i;
       }

    // this is out of bound, the compiler could complain or segmentation fault
    // or undermined behavior will likely follow at execution time
    //v[LENGTH]=1; 
    
    for(i=0; i<LENGTH; i++)
       {
       printf("%d ", v[i]);
       }
    printf("\n");

    times2(v);

    for(i=0; i<LENGTH; i++)
       {
       printf("%d ", v[i]);
       }
    printf("\n");

    //---------------- matrix

    for(i=0; i<LENGTH; i++)
       {
       for(j=0; j<LENGTH; j++)
          {
          matrix[i][j]=i+j;
          }
       }

    times2matrix(matrix);

    //-------------  strings

    // string initialized by using char
    name[0]='c';
    name[1]='i';
    name[2]='a';
    name[3]='o';
    for(i=0; i<20; i++)
       {
       printf("%c ", name[i]);
       }
    printf("\n");
    printf("%s\n", name);

    // string initialized as a string
    strcpy(name, "riciao");
    for(i=0; i<20; i++)
       {
       printf("%c ", name[i]);
       }
    printf("\n");
    printf("%s\n", name);

    return EXIT_SUCCESS;
    }


