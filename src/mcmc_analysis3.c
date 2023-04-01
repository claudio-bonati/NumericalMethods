#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/random.h"

#define STRING_LENGTH 50

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


// compute the jacknife samples of the Binder cumulant U
void binderUjack(double *datajack, double const * const data, long int numberofbins, int binsize)
  {
  long int i;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double x2tot, x4tot, x2, x4;

  x2tot=0.0;
  x4tot=0.0;
  for(i=0; i<sampleeff; i++)
     {
     x2tot+=pow(data[i],2.0);
     x4tot+=pow(data[i],4.0);
     }

  for(i=0; i<numberofbins; i++)
     {
     x2=x2tot;
     x4=x4tot;

     for(j=0; j<binsize; j++)
        {
        x2-=pow(data[i*binsize+j],2.0);
        x4-=pow(data[i*binsize+j],4.0);
        }
   
     x2/=(double)((numberofbins-1)*binsize);
     x4/=(double)((numberofbins-1)*binsize);
    
     datajack[i]=x4/(x2*x2);
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, ris, err;
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
    datajack=(double *)malloc((unsigned long int)numberofbins*sizeof(double));
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata(datafile, therm, sampleeff, data);

    // compute jackknife resamplings
    binderUjack(datajack, data, numberofbins, binsize);

    // compute average
    ris=0.0;
    for(i=0; i<numberofbins; i++)
       {
       ris+=datajack[i];
       }
    ris/=(double)numberofbins;

    // compute error
    err=0.0;
    for(i=0; i<numberofbins; i++)
       {
       err+=pow(ris-datajack[i], 2.0);
       }
    // this corrects for a factor that is irrelevant but we leave it just for clarity
    err*=(double)(numberofbins-1);
    err/=(double)numberofbins;
    err=sqrt(err);

    // free data arrays
    free(data);
    free(datajack);

    printf("%f %f\n", ris, err);

    return EXIT_SUCCESS;
    }

