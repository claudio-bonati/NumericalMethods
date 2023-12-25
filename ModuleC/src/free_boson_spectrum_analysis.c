#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize,
                 int Nt)
  {
  long int i, r, t, numcol;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double c_tot, cp1_tot;
  double c, cp1;

  // number of columns of datafile
  numcol=(Nt/4);

  for(r=0; r<numberofbins*(numcol-1); r++)
     {
     datajack[r]=0.0;
     }

  for(t=0; t<(Nt/4)-1; t++)
     {
     c_tot=0.0;
     cp1_tot=0.0;

     for(i=0; i<sampleeff; i++)
        {
        c_tot+=data[numcol*i+t];
        cp1_tot+=data[numcol*i+t+1];
        }

     for(i=0; i<numberofbins; i++)
        {
        c=c_tot;
        cp1=cp1_tot;

        for(j=0; j<binsize; j++)
           {
           r=i*binsize+j;
           c-=data[numcol*r+t];
           cp1-=data[numcol*r+t+1];
           }

        c/=(double)((numberofbins-1)*binsize);
        cp1/=(double)((numberofbins-1)*binsize);
  
        datajack[(numcol-1)*i + t]=log(c/cp1);
        }
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j, Nt, numcol;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, *ris, *err;
    double hatm;
    char datafile[STRING_LENGTH];

    if(argc != 6)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize Nt hatm datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  Nt = number of temporal steps\n");
      fprintf(stdout, "  hatm = a*m\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  effective mass at various time distances\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);
      Nt=atoi(argv[3]);
      hatm=atof(argv[4]);

      if(strlen(argv[5]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        strcpy(datafile, argv[5]);
        }
      }

    if(Nt<=1)
      {
      fprintf(stderr, "'Nt' must be larger than 1\n");
      return EXIT_FAILURE;
      }

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // number of columns of the file
    numcol=(Nt/4); 

    // determine the length of the file
    sample=linecounter_mc(datafile, numcol);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(numcol*sampleeff)*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(numberofbins*(numcol-1))*sizeof(double)); // in the jacknife we do not need the last column
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate ris, err
    ris=(double *)malloc((unsigned long int)(numcol-1)*sizeof(double));
    if(ris==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    err=(double *)malloc((unsigned long int)(numcol-1)*sizeof(double));
    if(err==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, numcol);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize, Nt);

    // compute average
    for(j=0; j<numcol-1; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[(numcol-1)*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<numcol-1; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[(numcol-1)*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    for(j=0; j<numcol-1; j++)
       {
       printf("%.12f %.12f\n", ris[j]/hatm, err[j]/hatm);
       }

    // free data arrays
    free(data);
    free(datajack);
    free(ris);
    free(err);

    return EXIT_SUCCESS;
    }


