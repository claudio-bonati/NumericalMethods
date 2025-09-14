#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// main
int main(void)
    {
    int x, i;

    x=0;
    for(i=0; i<4; i++)  // i++ stands for i=i+1, 
                        // i.e. increment i by one at each iteration
                        // pay attention to the difference between
                        // i++ and ++i [irrelevant in the present context]
                        // if i=1 then 
                        // y=i++ gives y=1, i=2
                        // y=++i dives y=2, i=2
       {
       x=x+2;  // we can also write x+=2
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

    printf("\n");

//--------------------------

    x=3;
    if(x<2)
      {
      printf("x<2\n");
      }
    else if(x==2) // pay attention to the "double" equal, 
           {      // with a single = the result would be complely different
                  // sometimes Yoda writing is suggested to avoid errors: 2==x
                  // will complain if used with a single =
           printf("x=2\n");
           }
         else if(x>2 && x<5)
                {
                printf("2<x<5\n");
                }
              else
                {
                printf("x>5\n");
                } 

    printf("\n");

//--------------------------

    x=0;
    switch(x)  // in general switch(expression)
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


