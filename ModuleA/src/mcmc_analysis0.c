#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// main
int main(int argc, char **argv)
    {
    int therm, j, maxt;
    long int i, sample, sampleeff;
    char datafile[STRING_LENGTH];
    double *data, corr, corr0, average;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm maxt datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  maxt = maximum time difference for \n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed (single column!)\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  t   <x(i)x(i+t)> for t=0 up to maxt\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      maxt=atoi(argv[2]);

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

    if(maxt<=0)
      {
      fprintf(stderr, "'maxt' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_sc(datafile);

    // effective number of data
    sampleeff=sample-therm;

    // allocate data array
    data=(double *)malloc((unsigned long int)sampleeff*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_sc(datafile, therm, sampleeff, data);

    // data average
    average=0.0;
    for(i=0; i<sampleeff; i++)
       {
       average+=data[i];
       }
    average/=(double)sampleeff;

    // compute and print correlators
    corr0=0.0;
    for(i=0; i<sampleeff; i++)
       {
       corr0+=(data[i]-average)*(data[i]-average);
       }
    corr0/=(double)sampleeff;
    printf("%d  %f\n", 0, 1.0);

    for(j=1; j<maxt; j++)
       {
       corr=0.0;

       for(i=0; i<sampleeff-j; i++)
          {
          corr+=(data[i]-average)*(data[i+j]-average);
          }

       corr/=(double)(sampleeff-j);
       printf("%d  %f\n", j, corr/corr0);
       }

    // free data array
    free(data);

    return EXIT_SUCCESS;
    }

