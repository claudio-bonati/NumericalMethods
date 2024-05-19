#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define LENGTH 10  // macro 
// if this number is too large one gest in execution
// Segmentation fault (core dumped)
// and dynamic memory allocation is required


void times2(int *x) // also "void times2(int x[LENGTH])" is ok
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
// and in fact it would be better: some compilers would detect the problem in
// the following function and complain
void times2wrong(int *x)
   {
   int i;
 
   for(i=0; i<LENGTH+10; i++)
      {
      x[i]*=2;
      }
   }


void times2matrix(int x[LENGTH][LENGTH])
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
    int vv[LENGTH][LENGTH];
    char name[20];

    for(i=0; i<LENGTH; i++)
       {
       v[i]=i;
       }

    // this is out of bound the compiler will complain or segmentation fault
    // will likely follows at execution time
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

    ///////

    for(i=0; i<LENGTH; i++)
       {
       for(j=0; j<LENGTH; j++)
          {
          vv[i][j]=i+j;
          }
       }

    times2matrix(vv);

    ////////

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


