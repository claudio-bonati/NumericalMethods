#include<stdio.h>
#include<stdlib.h>

// main
int main(int argc, char **argv)
    {
    int a, b, ris;

    if(argc != 3)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s a b\n\n", argv[0]);
      fprintf(stdout, "  a = integer\n");
      fprintf(stdout, "  b = integer\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  a*b\n\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      a=atoi(argv[1]);
      b=atoi(argv[2]);
      }


    ris=a*b;

    printf("%d*%d=%d\n", a, b, ris);

    return EXIT_SUCCESS;
    }


