#include<stdio.h>
#include<stdlib.h>

#include"../include/read_data.h"

//-------  SINGLE COLUMN ------- 

// determine the length of the file (single column!)
long int linecounter_sc(char const * const filename)
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


// initialize data (single column!)
void readdata_sc(char const * const filename, int therm, long int sampleeff, double *data)
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

