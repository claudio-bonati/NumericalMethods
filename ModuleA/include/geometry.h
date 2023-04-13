#ifndef GEOMETRY_H

// cartesian coordinates -> lexicographic index
// lattice L^{dim}
void cart_to_lex(long int * restrict lex, int const * const restrict cartcoord, int L, int dim);

// lexicographic index -> cartesian coordinates
// lattice L^{dim}
void lex_to_cart(int * restrict cartcoord, long int lex, int L, int dim);

// initialize geometry
// nnp[i*volume+r]= next neighbor in positive "i" direction of site r 
// nnm[i*volume+r]= next neighbor in negative "i" direction of site r 
void init_neighbors(long int * restrict nnp, long int * restrict nnm, int L, int dim);

void test_geometry(long int const * const nnp, long int const * const nnm, int L, int dim);


#endif
