#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/geometry.h"

// cartesian coordinates -> lexicographic index
// lattice L^{dim}
void cart_to_lex(long int * restrict lex, int const * const restrict cartcoord, int L, int dim)
  {
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<dim; i++)
     {
     ris+=cartcoord[i]*aux;
     aux*=L;  
     }

  // ris = cartcoord[0]
  //      +cartcoord[1]*L
  //      +cartcoord[2]*L*L
  //      +...
  //      +cartcoord[dim-1]*L^{dim-1}

  *lex=ris;
  }


// lexicographic index -> cartesian coordinates
// lattice L^{dim}
void lex_to_cart(int * restrict cartcoord, long int lex, int L, int dim)
  {
  int i;
  long aux;

  aux=1;
  for(i=0; i<dim-1; i++)
     {
     aux*=L;
     }
  // aux=pow(L,dim-1) but pow uses float

  for(i=dim-1; i>=0; i--)
     {
     cartcoord[i]=(int) (lex/aux);

     lex-=aux*cartcoord[i];
     aux/=L;  
     }
  }


// inline function defined in include/geometry.h
//
// nnp[dirgeo(r, i, volume)] is the neighbor of "r" in positive direction "i" on a lattice of volume "volume"
// nnm[dirgeo(r, i, volume)] is the neighbor of "r" in negative direction "i" on a lattice of volume "volume"
long int dirgeo(long int lex, int i, long int volume);


// initialize geometry
// nnp[dirgeo(r,i,volume)]= next neighbor in positive "i" direction of site r 
// nnm[dirgeo(r,i,volume)]= next neighbor in negative "i" direction of site r 
void init_neighbors(long int * restrict nnp, 
                    long int * restrict nnm, 
                    int L, 
                    int dim)
  {
  int i, value, valuep, valuem;
  long r, rm, rp, volume;
  int *cartcoord;

  // allocate caresian coordinate vector
  cartcoord=(int *)malloc((unsigned long int)dim*sizeof(int));
  if(cartcoord==NULL)
    {
    fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // pow is a float function, so we compute the power by hand
  volume=1;
  for(i=0; i<dim; i++)
     {
     volume*=L; 
     }

  // initialize next neighbors
  for(r=0; r<volume; r++)
     {
     lex_to_cart(cartcoord, r, L, dim);

     for(i=0; i<dim; i++)
        {
        value=cartcoord[i];

        valuep=value+1;
        if(valuep >= L)
          {
          valuep-=L;
          }
        cartcoord[i]=valuep;
        cart_to_lex(&rp, cartcoord, L, dim);
        nnp[dirgeo(r, i, volume)]=rp;

        valuem=value-1;
        if(valuem<0)
          {
          valuem+=L;
          }
        cartcoord[i]=valuem;
        cart_to_lex(&rm, cartcoord, L, dim);
        nnm[dirgeo(r, i, volume)]=rm;

        cartcoord[i]=value;
        }
     } // end of loop on r

  // free caresian coordinate vector
  free(cartcoord);

  #ifdef DEBUG
    test_geometry(nnp, nnm, L, dim);
  #endif
  }  


void test_geometry(long int const * const restrict nnp, 
                   long int const * const restrict nnm, 
                   int L, 
                   int dim)
  {
  int i, dir, dir1, *cartcoord;
  long r, volume, r_test, r_test1, r_test2;

  // allocate caresian coordinate vector
  cartcoord=(int *)malloc((unsigned long int)dim*sizeof(int));
  if(cartcoord==NULL)
    {
    fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // pow is a float function, so we compute the power by hand
  volume=1;
  for(i=0; i<dim; i++)
     {
     volume*=L; 
     }

  // test of lex <-> cart
  for(r=0; r < volume; r++)
     {
     lex_to_cart(cartcoord, r, L, dim);
     cart_to_lex(&r_test, cartcoord, L, dim);

     if(r != r_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of nnp <-> nnm
  for(r=0; r < volume; r++)
     {
     for(dir=0; dir<dim; dir++)
        {
        r_test=nnp[dirgeo(r, dir, volume)];
        r_test1=nnm[dirgeo(r_test, dir, volume)];

        if(r != r_test1)
          {
          fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

  // test of nnp nnp <-> nnm nnm along different directions
  for(r=0; r < volume; r++)
     {
     for(dir=0; dir<dim; dir++)
        {
        r_test1=nnp[dirgeo(r, dir, volume)];

        for(dir1=0; dir1<dim; dir1++)
           {
           r_test2=nnp[dirgeo(r_test1, dir1, volume)];

           r_test=nnm[dirgeo(r_test2, dir, volume)];
           r_test2=nnm[dirgeo(r_test, dir1, volume)];

           if(r != r_test2)
             {
             fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }
           }
        }
     }

  fprintf(stdout, "Geometry test passed! (%s, %d)\n", __FILE__, __LINE__);

  free(cartcoord);
  }

