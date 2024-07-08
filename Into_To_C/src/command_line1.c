#include<stdio.h>
#include<stdlib.h>

// main
int main(int argc, char **argv) //<---- in previous examples it was int main(void)
    {
    int a, b, ris;

    if(argc != 3) // 3 means that two arguments are required
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s a b\n\n", argv[0]); // argv[0] is the executable name
      fprintf(stdout, "  a = integer\n");
      fprintf(stdout, "  b = integer\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  a*b\n\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      a=atoi(argv[1]); // atoi=library function to convert string to integer
      b=atoi(argv[2]);
      }


    ris=a*b;

    printf("%d*%d=%d\n", a, b, ris);

    return EXIT_SUCCESS;
    }


