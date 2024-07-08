#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define STRING_LENGTH 50

// main
int main(int argc, char **argv)
    {
    double a, b, ris;
    char datafile[STRING_LENGTH];
    FILE *fp;

    if(argc != 4) // 4 -> 3 arguments are expected
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s a b filename\n\n", argv[0]);
      fprintf(stdout, "  a = real number\n");
      fprintf(stdout, "  b = real number\n");
      fprintf(stdout, "  filename = name of the file on which the result is written\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  sqrt(a*b)\n\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      a=atof(argv[1]); // atof = string to double 
      b=atof(argv[2]);

      // filename could be too long, we have to check to prevent errors
      if(strlen(argv[3]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        strcpy(datafile, argv[3]); // copy argv[3] in the string datafile
        }
      }

    // compute the result
    ris=sqrt(a*b);

    // open data file for writing 
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // print the result on the file in different ways
    fprintf(fp, "sqrt(%lf*%lf)=%lf\n", a, b, ris);
    fprintf(fp, "with more places: sqrt(%.8lf*%.8lf)=%.12lf\n", a, b, ris);
    fprintf(fp, "with too many places: sqrt(%.25lf*%.25lf)=%.25lf\n", a, b, ris);
    fprintf(fp, "check the last results on https://www.wolframalpha.com\n");

    // close the file
    fclose(fp);

    return EXIT_SUCCESS;
    }


