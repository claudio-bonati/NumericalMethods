#include<stdio.h>
#include<stdlib.h>
#include<string.h>


// main
int main(void)
    {
    int i, err;
    double x[10];
    int xi[10];
    char datafile[20]; // file name
    FILE *fp; // pointer to file

    strcpy(datafile, "data.dat"); // file name initialized with a string

    // open data file for writing (overwrite the previous content)
    // if "a" is used instead of "w" data are appended, without deleting
    // the previous content  
    fp=fopen(datafile, "w");
    //// it is safer (although not mandatory) to check that 
    //// everything worked properly, as done below
    //if(fp==NULL)
    //  {
    //  fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    //  return EXIT_FAILURE;
    //  }

    for(i=0; i<10; i++)
       {
       fprintf(fp, "%lf  %d\n", ((double) i)/10.0, i);
       //// we are neglecting the possibility of errors in frprintf. It would be
       //// more correct to use the form below
       //err = fprintf(fp, "%lf  %d\n", ((double) i)/10.0, i);
       //if(err!=2)
       //  {
       //  fprintf(stderr, "Error in printf (%s, %d)\n", __FILE__, __LINE__);
       //  exit(EXIT_FAILURE);
       //  }
       }

    fclose(fp); // close file 

    // open data file for reading
    fp=fopen(datafile, "r");
    //// it is safer to check that everything worked properly, as done below
    //if(fp==NULL)
    //  {
    //  fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    //  return EXIT_FAILURE;
    //  }
 
    for(i=0; i<10; i++)
       {
       err=fscanf(fp, "%lf %d", &(x[i]), &(xi[i])); // neglecting "err" would result in a warning
                                                    // if everything is ok "err" should be equal to 2
                                                    // (in general fscanf returns the number of objects read)
       if(err!=2)
         {
         fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__); // print on the standard error the file 
                                                                           // and line at which the error occurred
         exit(EXIT_FAILURE); // this abort the main, returning an error code
         }
       }

    fclose(fp); // close file

    for(i=0; i<10; i++)
       {
       printf("%lf ", x[i]);
       }
    printf("\n");

    for(i=0; i<10; i++)
       {
       printf("%d ", xi[i]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }


