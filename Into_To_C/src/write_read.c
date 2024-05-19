#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


// main
int main(void)
    {
    int i, err;
    double x[10];
    char datafile[20];
    FILE *fp; // pointer to file

    strcpy(datafile, "data.dat");

    // open data file for writing (overwrite the previous content)
    fp=fopen(datafile, "w");
  
    for(i=0; i<10; i++)
       {
       fprintf(fp, "%lf\n", ((double) i)/10.0);
       // we are neglecting the possibility of errors in frprintf. It would be more correct to use err=fprintf...
       }

    fclose(fp);

    // open data file for reading
    fp=fopen(datafile, "r");
  
    for(i=0; i<10; i++)
       {
       err=fscanf(fp, "%lf\n", &(x[i]));
       if(err!=1)
         {
         fprintf(stderr, "Error in scanf (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);

    for(i=0; i<10; i++)
       {
       printf("%lf ", x[i]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }


