#include<stdio.h>
#include<stdlib.h>

// main
int main()
    {
    int *vec;
    long int r, length=1000;

    // allocate the vector. If length is too big this will fail in execution (try it!)
    vec=(int *)malloc((unsigned long int)(length)*sizeof(int));
    if(vec == NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialization with +1 and -1
    for(r=0; r<length; r++)
       {
       if(r%2==0)
         {
         vec[r]=1;
         }
       else
         {
         vec[r]=-1;
         }
       }

    printf("%d\n", vec[126]);

    // free the memory
    free(vec);

    // if not commented this line typically produces "Segmentation fault (core dumped)" since we are accessing a region of memory that is no more legitimate (try it!)
    // printf("%d\n", vec[126]);

    return EXIT_SUCCESS;
    }


