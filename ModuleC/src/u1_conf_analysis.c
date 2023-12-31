#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/read_data.h"

#define MIN(a,b) (((a)<(b))?(a):(b))  // min of two numbers

#define STRING_LENGTH 50

void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize,
                 int Nt, 
                 int Ns)
  {
  long int i, s, t, r, aux, numcol;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double c_tot;
  double c;

  // number of columns of datafile
  numcol=(Ns/4)*MIN((Nt/4),8);

  for(r=0; r<numberofbins*numcol; r++)
     {
     datajack[r]=0.0;
     }

  for(s=0; s<(Ns/4); s++)
     {
     for(t=0; t<MIN((Nt/4),8); t++)
        {
        aux=s*MIN((Nt/4),8)+t;

        c_tot=0.0;

        for(i=0; i<sampleeff; i++)
           {
           c_tot+=data[numcol*i+aux];
           }

        for(i=0; i<numberofbins; i++)
           {
           c=c_tot;

           for(j=0; j<binsize; j++)
              {
              r=i*binsize+j;
              c-=data[numcol*r+aux];
              }

           c/=(double)((numberofbins-1)*binsize);
  
           datajack[numcol*i + aux] = -log(c)/(double)(t+1);
           }
        }
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j, k, Nt, Ns, numcol;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, *ris, *err;
    char datafile[STRING_LENGTH];

    if(argc != 6)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize Nt Ns datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  Nt = number of temporal sites\n");
      fprintf(stdout, "  Ns = number of spatial sites\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  effective potential computed from Wilson loops of\n");
      fprintf(stdout, "  different time extent\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);
      Nt=atoi(argv[3]);
      Ns=atoi(argv[4]);

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

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // number of columns of the file
    numcol=MIN((Nt/4),8)*(Ns/4); 

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
    datajack=(double *)malloc((unsigned long int)(numberofbins*numcol)*sizeof(double)); 
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate ris, err
    ris=(double *)malloc((unsigned long int)(numcol)*sizeof(double));
    if(ris==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    err=(double *)malloc((unsigned long int)(numcol)*sizeof(double));
    if(err==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, numcol);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize, Nt, Ns);

    // compute average
    for(j=0; j<numcol; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[numcol*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<numcol; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[numcol*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    for(j=0; j<Ns/4; j++)
       {
       printf("%d ",j+1);
       for(k=0; k<MIN(Nt/4,8); k++)
          {
          printf("%.12f %.12f ", ris[j*MIN((Nt/4),8)+k], err[j*MIN((Nt/4),8)+k]);
          }
       printf("\n");
       }

    // free data arrays
    free(data);
    free(datajack);
    free(ris);
    free(err);

    return EXIT_SUCCESS;
    }


