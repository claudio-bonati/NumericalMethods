#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// main
int main(void)
    {
    int x, i;

    x=0;
    for(i=0; i<4; i++)
       {
       x=x+2;
       }
    printf("1)  %d\n", x);

    i=0;
    x=0;
    while(i<4)
      {
      x=x+2;
      i++;
      }
    printf("2)  %d\n", x);

    i=0;
    x=0;
    do
      {
      x=x+2;
      i++;
      }
    while(i<4);

    printf("3)  %d\n", x);

/////
    printf("\n");

    x=3;
    if(x<=2)
      {
      printf("x<=2\n");
      }
    else if(x>2 && x<5)
           {
           printf("2<x<5\n");
           }
         else
           {
           printf("x>5\n");
           } 

/////
    printf("\n");


    x=3;
    switch(x)
      {
      case 0: printf("x=0\n");
              break;
      case 1: printf("x=1\n");
              break;
      case 2: printf("x=2\n");
              break;
      case 3: printf("x=3\n");
              break;
      default: printf("something else\n");
      }

    return EXIT_SUCCESS;
    }


