#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/random.h"

#define STRING_LENGTH 50
#define MAXBOOT 100

// determine the length of the file
long int linecounter(char filename[STRING_LENGTH])
  {
  int err;
  long int sample;
  double tmp;
  FILE *fp;

  // open data file
  fp=fopen(filename, "r");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  
  // count lines of datafile
  err=1;
  sample=0;
  while(err==1)
    {
    sample++;
    err=fscanf(fp, "%lf", &tmp);
    } 
  sample--; // the last one has to be removed since err!=1

  // close datafile
  fclose(fp);

  return sample;
  }


// initialize data
void readdata(char filename[STRING_LENGTH], int therm, long int sampleeff, double *data)
  {
  int err;
  long int i;
  double tmp;
  FILE *fp;

  // open data file
  fp=fopen(filename, "r");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // read thermalization
  for(i=0; i<therm; i++)
     {
     err=fscanf(fp, "%lf", &tmp);
     if(err!=1)
       {
       fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // intialize data array
  for(i=0; i<sampleeff; i++)
     {
     err=fscanf(fp, "%lf", &tmp);
     if(err!=1)
       {
       fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     data[i]=tmp;
     }

  fclose(fp);
  }


// resample
void resample(double *boot, const double * const data, long int numberofbins, int binsize)
  {
  long int i, ib;
  double tmp;
  int j;

  for(i=0; i<numberofbins*binsize; i++)
     {
     boot[i]=data[i];
     }

  for(i=0; i<numberofbins; i++)
     {
     tmp=(double)numberofbins*myrand();
     ib=(long int) tmp;
     for(j=0; j<binsize; j++)
        {
        boot[i*binsize+j]=data[ib*binsize+j];
        }
     }
  }


// compute the Binder cumulant U
double binderU(double *boot, long int sampleeff)
  {
  long int i;
  double x2, x4;

  x2=0.0;
  x4=0.0;

  for(i=0; i<sampleeff; i++)
     {
     x2+=pow(boot[i],2.0);
     x4+=pow(boot[i],4.0);
     }

  x2/=(double)sampleeff;
  x4/=(double)sampleeff;

  return x4/(x2*x2);
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, boot;
    long int sample, numberofbins, sampleeff;
    double *data, *databoot, U[MAXBOOT], ris, err;
    char datafile[STRING_LENGTH];

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed (single column!)\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <x^4>/<x^2>^2 and it error, computed using binning and bootstrap\n");

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
    sample=linecounter(datafile);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)sampleeff*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    databoot=(double *)malloc((unsigned long int)sampleeff*sizeof(double));
    if(databoot==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata(datafile, therm, sampleeff, data);

    for(boot=0; boot<MAXBOOT; boot++)
       {
       // resample
       resample(databoot, data, numberofbins, binsize);

       U[boot]=binderU(databoot, sampleeff);
       }

    // compute average and its std
    ris=0.0;
    for(boot=0; boot<MAXBOOT; boot++)
       {
       ris+=U[boot];
       }
    ris/=(double)MAXBOOT;

    err=0.0;
    for(boot=0; boot<MAXBOOT; boot++)
       {
       err+=pow((ris-U[boot]),2.0);
       }
    err/=(double)MAXBOOT;
    err=sqrt(err);
  

    // free data arrays
    free(data);
    free(databoot);

    printf("%f %f\n", ris, err);

    return EXIT_SUCCESS;
    }

