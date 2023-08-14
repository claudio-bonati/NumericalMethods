#ifndef GEOMETRY_H
#define GEOMETRY_H

// cartesian coordinates -> lexicographic index
// lattice L^{dim}
void cart_to_lex(long int * restrict lex, int const * const restrict cartcoord, int L, int dim);

// lexicographic index -> cartesian coordinates
// lattice L^{dim}
void lex_to_cart(int * restrict cartcoord, long int lex, int L, int dim);

// nnp[dirgeo(r, i, volume)] is the neighbor of "r" in positive direction "i" on a lattice of volume "volume"
// nnm[dirgeo(r, i, volume)] is the neighbor of "r" in negative direction "i" on a lattice of volume "volume"
inline long int dirgeo(long int lex, int i, long int volume)
  {
  return i*volume+lex;
  }

// initialize geometry
// nnp[dirgeo(r,i,volume)]= next neighbor in positive "i" direction of site r 
// nnm[dirgeo(r,i,volume)]= next neighbor in negative "i" direction of site r 
void init_neighbors(long int * restrict nnp, long int * restrict nnm, int L, int dim);

void test_geometry(long int const * const restrict nnp, 
                   long int const * const restrict nnm, 
                   int L, 
                   int dim);


#endif
