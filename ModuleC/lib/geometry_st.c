#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/geometry_st.h"

//#define DEBUG

// cartesian coordinates -> lexicographic index
// lattice Nt*Ns^{stdim-1}
void cart_to_lex_st(long int * restrict lex, int const * const restrict cartcoord, int Nt, int Ns, int stdim)
  {
  int i;
  long ris, aux;

  ris=0;
  aux=1;

  ris+=cartcoord[0]*aux;
  aux*=Nt;  
  for(i=1; i<stdim; i++)
     {
     ris+=cartcoord[i]*aux;
     aux*=Ns;  
     }

  // ris = cartcoord[0]
  //      +cartcoord[1]*Nt
  //      +cartcoord[2]*Nt*Ns
  //      +...
  //      +cartcoord[stdim-1]*Nt*Ns^{stdim-2}

  *lex=ris;
  }


// lexicographic index -> cartesian coordinates
// lattice Nt*Ns^{stdim-1}
void lex_to_cart_st(int * restrict cartcoord, long int lex, int Nt, int Ns, int stdim)
  {
  int i;
  long aux;

  aux=Nt;
  for(i=1; i<stdim-1; i++)
     {
     aux*=Ns;
     }
  // aux=Nt*pow(Ns,stdim-2) but pow uses float

  for(i=stdim-1; i>=2; i--)
     {
     cartcoord[i]=(int) (lex/aux);

     lex-=aux*cartcoord[i];
     aux/=Ns;  
     }

  cartcoord[1]=(int) (lex/aux);
  lex-=aux*cartcoord[i];
  aux/=Nt;  

  cartcoord[0]=(int) (lex/aux);
  }


// inline function defined in include/geometry.h
//
// nnp[dirgeo(r, i, volume)] is the neighbor of "r" in positive direction "i" on a lattice of volume "volume"
// nnm[dirgeo(r, i, volume)] is the neighbor of "r" in negative direction "i" on a lattice of volume "volume"
long int dirgeo(long int lex, int i, long int volume);


// initialize geometry
// nnp[dirgeo(r,i,volume)]= next neighbor in positive "i" direction of site r 
// nnm[dirgeo(r,i,volume)]= next neighbor in negative "i" direction of site r 
void init_neighbors_st(long int * restrict nnp, 
                       long int * restrict nnm, 
                       int Nt,
                       int Ns,  
                       int stdim)
  {
  int i, value, valuep, valuem;
  long r, rm, rp, stvolume;
  int *cartcoord, *size;

  // allocate caresian coordinate vector
  cartcoord=(int *)malloc((unsigned long int)stdim*sizeof(int));
  if(cartcoord==NULL)
    {
    fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // lattice sizes
  size=(int *)malloc((unsigned long int)stdim*sizeof(int));
  if(size==NULL)
    {
    fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  size[0]=Nt;
  for(i=1; i<stdim; i++)
     {
     size[i]=Ns;
     }

  // pow is a float function, so we compute the power by hand
  stvolume=Nt;
  for(i=1; i<stdim; i++)
     {
     stvolume*=Ns; 
     }

  // initialize next neighbors
  for(r=0; r<stvolume; r++)
     {
     lex_to_cart_st(cartcoord, r, Nt, Ns, stdim);

     for(i=0; i<stdim; i++)
        {
        value=cartcoord[i];

        valuep=value+1;
        if(valuep >= size[i])
          {
          valuep-=size[i];
          }
        cartcoord[i]=valuep;
        cart_to_lex_st(&rp, cartcoord, Nt, Ns, stdim);
        nnp[dirgeo(r, i, stvolume)]=rp;

        valuem=value-1;
        if(valuem<0)
          {
          valuem+=size[i];
          }
        cartcoord[i]=valuem;
        cart_to_lex_st(&rm, cartcoord, Nt, Ns, stdim);
        nnm[dirgeo(r, i, stvolume)]=rm;

        cartcoord[i]=value;
        }
     } // end of loop on r

  // free caresian coordinate vector
  free(cartcoord);
  free(size);

  #ifdef DEBUG
    test_geometry_st(nnp, nnm, Nt, Ns, stdim);
  #endif
  }  


void test_geometry_st(long int const * const restrict nnp, 
                      long int const * const restrict nnm, 
                      int Nt,
                      int Ns, 
                      int stdim)
  {
  int i, dir, dir1, *cartcoord;
  long r, stvolume, r_test, r_test1, r_test2;

  // allocate caresian coordinate vector
  cartcoord=(int *)malloc((unsigned long int)stdim*sizeof(int));
  if(cartcoord==NULL)
    {
    fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // pow is a float function, so we compute the power by hand
  stvolume=Nt;
  for(i=1; i<stdim; i++)
     {
     stvolume*=Ns; 
     }

  // test of lex <-> cart
  for(r=0; r < stvolume; r++)
     {
     lex_to_cart_st(cartcoord, r, Nt, Ns, stdim);
     cart_to_lex_st(&r_test, cartcoord, Nt, Ns, stdim);

     if(r != r_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of nnp <-> nnm
  for(r=0; r < stvolume; r++)
     {
     for(dir=0; dir<stdim; dir++)
        {
        r_test=nnp[dirgeo(r, dir, stvolume)];
        r_test1=nnm[dirgeo(r_test, dir, stvolume)];

        if(r != r_test1)
          {
          fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

  // test of nnp nnp <-> nnm nnm along different directions
  for(r=0; r < stvolume; r++)
     {
     for(dir=0; dir<stdim; dir++)
        {
        r_test1=nnp[dirgeo(r, dir, stvolume)];

        for(dir1=0; dir1<stdim; dir1++)
           {
           r_test2=nnp[dirgeo(r_test1, dir1, stvolume)];

           r_test=nnm[dirgeo(r_test2, dir, stvolume)];
           r_test2=nnm[dirgeo(r_test, dir1, stvolume)];

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

