#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

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
  double c_tot, cp1_tot, x2_tot;
  double c, cp1, x2;

  // number of columns of datafile
  numcol=4*(Nt/4)+1;

  for(r=0; r<numberofbins*(numcol-5); r++)
     {
     datajack[r]=0.0;
     }

  // x x correlator
  for(t=0; t<(Nt/4)-1; t++)
     {
     c_tot=0.0;
     cp1_tot=0.0;

     for(i=0; i<sampleeff; i++)
        {
        c_tot+=data[numcol*i+1+4*t];
        cp1_tot+=data[numcol*i+1+4*(t+1)];
        }

     for(i=0; i<numberofbins; i++)
        {
        c=c_tot;
        cp1=cp1_tot;

        for(j=0; j<binsize; j++)
           {
           r=i*binsize+j;
           c-=data[numcol*r+1+4*t];
           cp1-=data[numcol*r+1+4*(t+1)];
           }

        c/=(double)((numberofbins-1)*binsize);
        cp1/=(double)((numberofbins-1)*binsize);
  
        datajack[(numcol-5)*i + 4*t + 0]=log(c/cp1);
        }
     }

  // x2 x2 correlator
  for(t=0; t<(Nt/4)-1; t++)
     {
     c_tot=0.0;
     cp1_tot=0.0;
     x2_tot=0.0;

     for(i=0; i<sampleeff; i++)
        {
        c_tot+=data[numcol*i+2+4*t];
        cp1_tot+=data[numcol*i+2+4*(t+1)];
        x2_tot+=data[numcol*i+0];
        }

     for(i=0; i<numberofbins; i++)
        {
        c=c_tot;
        cp1=cp1_tot;
        x2=x2_tot;

        for(j=0; j<binsize; j++)
           {
           r=i*binsize+j;
           c-=data[numcol*r+2+4*t];
           cp1-=data[numcol*r+2+4*(t+1)];
           x2-=data[numcol*r+0];
           }

        c/=(double)((numberofbins-1)*binsize);
        cp1/=(double)((numberofbins-1)*binsize);
        x2/=(double)((numberofbins-1)*binsize);
 
        datajack[(numcol-5)*i + 4*t + 1]=log((c-x2*x2)/(cp1-x2*x2));
        }
     }

  // x3 x3 correlator
  for(t=0; t<(Nt/4)-1; t++)
     {
     c_tot=0.0;
     cp1_tot=0.0;

     for(i=0; i<sampleeff; i++)
        {
        c_tot+=data[numcol*i+3+4*t];
        cp1_tot+=data[numcol*i+3+4*(t+1)];
        }

     for(i=0; i<numberofbins; i++)
        {
        c=c_tot;
        cp1=cp1_tot;

        for(j=0; j<binsize; j++)
           {
           r=i*binsize+j;
           c-=data[numcol*r+3+4*t];
           cp1-=data[numcol*r+3+4*(t+1)];
           }

        c/=(double)((numberofbins-1)*binsize);
        cp1/=(double)((numberofbins-1)*binsize);
 
        datajack[(numcol-5)*i + 4*t + 2]=log(c/cp1);
        }
     }

  // A A correlator
  for(t=0; t<(Nt/4)-1; t++)
     {
     c_tot=0.0;
     cp1_tot=0.0;

     for(i=0; i<sampleeff; i++)
        {
        c_tot+=data[numcol*i+4+4*t];
        cp1_tot+=data[numcol*i+4+4*(t+1)];
        }

     for(i=0; i<numberofbins; i++)
        {
        c=c_tot;
        cp1=cp1_tot;

        for(j=0; j<binsize; j++)
           {
           r=i*binsize+j;
           c-=data[numcol*r+4+4*t];
           cp1-=data[numcol*r+4+4*(t+1)];
           }

        c/=(double)((numberofbins-1)*binsize);
        cp1/=(double)((numberofbins-1)*binsize);
 
        datajack[(numcol-5)*i + 4*t + 3]=log(c/cp1);
        }
     }

  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j, Nt, numcol;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, *ris, *err;
    double simbeta, eta;
    char datafile[STRING_LENGTH];

    if(argc != 6)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile simbeta Nt\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n");
      fprintf(stdout, "  simbeta = hbar*omega/(k_B T)\n");
      fprintf(stdout, "  Nt = number of temporal steps\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  effective mass computed from the correlators x-x, x^2-x^2, x^3-x^3, A-A (with A=x^3-(3/2)s)\n");
      fprintf(stdout, "  each line correspond to a different time, from zero to Nt/4-1\n\n");

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
    simbeta=atof(argv[4]);
    Nt=atoi(argv[5]);

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

    // eta
    eta=simbeta/Nt;

    // number of columns of the file
    numcol=4*(Nt/4)+1; 

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
    datajack=(double *)malloc((unsigned long int)(numberofbins*(numcol-5))*sizeof(double)); // in the jacknife we do not need the first and the last 4 columns
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate ris, err
    ris=(double *)malloc((unsigned long int)(numcol-5)*sizeof(double));
    if(ris==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    err=(double *)malloc((unsigned long int)(numcol-5)*sizeof(double));
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
    for(j=0; j<numcol-5; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[(numcol-5)*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<numcol-5; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[(numcol-5)*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    for(j=0; j<((Nt/4)-1); j++)
       {
       printf("%.12f %.12f ", ris[4*j]/eta,   err[4*j]/eta);
       printf("%.12f %.12f ", ris[4*j+1]/eta, err[4*j+1]/eta);
       printf("%.12f %.12f ", ris[4*j+2]/eta, err[4*j+2]/eta);
       printf("%.12f %.12f\n",ris[4*j+3]/eta, err[4*j+3]/eta);
       }

    // free data arrays
    free(data);
    free(datajack);
    free(ris);
    free(err);

    return EXIT_SUCCESS;
    }

