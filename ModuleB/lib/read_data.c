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
void readdata_sc(char const * const filename, int therm, long int sampleeff, double * restrict data)
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


// ---------------- MULTIPLE COLUMNS

// determine the length of the file with 'col' columns
long int linecounter_mc(char const * const filename, int col)
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

  if(sample % col !=0)
    {
    fprintf(stderr, "Something wrong in %s (%s, %d)\n", filename, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
   
  return sample/col;
  }


// initialize data from file with 'col' columns
// *data is structured as data[i*col+j]= j-th column of the i-th raw 
void readdata_mc(char const * const filename, int therm, long int sampleeff, double * restrict data, int col)
  {
  int err, j;
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
     for(j=0; j<col; j++)
        {
        err=fscanf(fp, "%lf", &tmp);
        if(err!=1)
          {
          fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

  // intialize data array
  for(i=0; i<sampleeff; i++)
     {
     for(j=0; j<col; j++)
        {
        err=fscanf(fp, "%lf", &tmp);
        if(err!=1)
          {
          fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }

        data[i*col+j]=tmp;
        }
     }

  fclose(fp);
  }


