#include<stdio.h>
#include<stdlib.h>

// main
int main(void)
    {
    int *vec;
    int **matrix;
    long int r, length=100;

    // allocate the vector. If length is too large this will fail in execution
    vec=(int *)malloc((unsigned long int)(length)*sizeof(int));
    // this check could be avoided but it is safer to used it
    if(vec == NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
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

    printf("%d\n", vec[length-1]);

    // if not commented this line typically produces 
    // "Segmentation fault (core dumped)" 
    // since we are accessing a region of memory that 
    // is not legitimate, but execution can also go on 
    // with unpredictable consequences.
    //
    // printf("%d\n", vec[length+10]);

    // free the memory of the vector
    free(vec);

    // if not commented this line typically produces 
    // "Segmentation fault (core dumped)" 
    // since we are accessing a region of memory that 
    // is NO MORE legitimate but execution can also go on 
    // with unpredictable consequences.
    //
    // printf("%d\n", vec[length-1]);

    //---------------------------

    // allocate the matrix
    matrix=(int **)malloc((unsigned long int)(length)*sizeof(int *));
    if(matrix == NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    else
      {
      for(r=0; r<length; r++)
         {
         matrix[r]=(int *)malloc((unsigned long int)(length)*sizeof(int));
         if(matrix[r] == NULL)
           {
           fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
           return EXIT_FAILURE;
           }
         }
      }

    // initialization of the first element
    matrix[0][0]=1;

    // free the memory of the matrix
    for(r=0; r<length; r++)
       {
       free(matrix[r]);
       }
    free(matrix);

    return EXIT_SUCCESS;
    }


