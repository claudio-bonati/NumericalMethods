#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/random.h"
#include"../include/read_data.h"

#define STRING_LENGTH 50

// main
int main(int argc, char **argv)
    {
    int j, therm, binsize, err;
    long int i, sample, numberofbins;
    double tmp, x, x2, xb2, x2b2, xtmp, x2tmp;
    char datafile[STRING_LENGTH];
    FILE *fp;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed (single column!)\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <x> and it error, <x^2> and its error\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);

      if(strlen(argv[3]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        strcpy(datafile, argv[3]);
        }
      }

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_sc(datafile);

    // initialize numberofbins
    numberofbins=(sample-therm)/binsize;

    // open data file [this is just a simple way to go back to the beginning of the file]
    fp=fopen(datafile, "r");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // read thermalization
    for(i=0; i<therm; i++)
       {
       err=fscanf(fp, "%lf", &tmp);
       if(err!=1)
         {
         fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
         return EXIT_FAILURE;
         }
       }

    // intialize average values
    x=0.0;
    x2=0.0;
    xb2=0.0;  // (binned x) squared
    x2b2=0.0; // (x^2 binned) squared
 
    for(i=0; i<numberofbins; i++)
       {
       xtmp=0.0;
       x2tmp=0.0;
 
       for(j=0; j<binsize; j++)
          {
          err=fscanf(fp, "%lf", &tmp);
          if(err!=1)
            {
            fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
            }
   
          xtmp+=tmp;
          x2tmp+=pow(tmp,2.0);
          }

       xtmp/=(double) binsize;
       x2tmp/=(double) binsize;

       x+=xtmp;
       xb2+=xtmp*xtmp;

       x2+=x2tmp;
       x2b2+=x2tmp*x2tmp;
       }

    x/=(double)numberofbins;
    xb2/=(double)numberofbins;
    x2/=(double)numberofbins;
    x2b2/=(double)numberofbins;

    printf("%f %f ", x, sqrt(xb2-x*x)/sqrt((double)numberofbins));
    printf("%f %f\n", x2, sqrt(x2b2-x2*x2)/sqrt((double)numberofbins));

    return EXIT_SUCCESS;
    }

