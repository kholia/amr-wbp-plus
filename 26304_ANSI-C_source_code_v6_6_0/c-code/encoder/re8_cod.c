#include "../include/amr_plus.h"
#include <assert.h> /* pre-release check */
#include <stdlib.h> /* abs() */
static void re8_compute_base_index( int x[], int *ka, long *I );
static void re8_compute_rank_of_permutation_and_sign_code( int *x, int *t, int *code );
/*--------------------------------------------------------------------------
  RE8_cod(x, n, I, k)
  MULTI-RATE INDEXING OF A POINT y in THE LATTICE RE8 (INDEX COMPUTATION)
  (i) x: point in RE8 (8-dimensional integer vector)
  (i) n: codebook number (*n is an integer defined in {0,2,3,4,..,n_max})
  (o) I: index of c (pointer to unsigned 16-bit word)
  (o) k: index of v (8-dimensional vector of binary indices) = Voronoi index
  note: the index I is defined as a 32-bit word, but only
  16 bits are required (long can be replaced by unsigned integer)
  --------------------------------------------------------------------------
 */
void RE8_cod( int x[], int *n, long *I, int k[] )
{
	int ka, c[8];
	/* decompose x as x = 2^r c + v, where r is an integer >=0, c is an element
     of Q0, Q2, Q3 or Q4, and v is an element of a Voronoi code in RE8
     (if r=0, x=c)
     this decomposition produces as a side-product the index k[] of v
     and the identifier ka of the absolute leader related to c
     the index of y is split into 2 parts :
     - the index I of c
     - the index k[] of v */
	RE8_vor( x, n, k, c, &ka );
	/* compute the index I (only if c is in Q2, Q3 or Q4) */
	if( *n > 0 ) {
		re8_compute_base_index( c, &ka, I );
	}
	return;
}
/*--------------------------------------------------------------
  re8_compute_base_index(x, ka, I)
  COMPUTE THE INDEX I of A LATTICE POINT x in  RE8
  (i) x : point in RE8 (8-dimensional integer vector)
  (i) ka: identifier of the absolute leader of x (scalar integer)
  (o) I : index of x (unsigned 16-bit word)
  note: the index I is defined as a 32-bit word, but only
  16 bits are required (long can be replaced by unsigned integer)
  ---------------------------------------------------------------
 */
void re8_compute_base_index( int x[], int *ka, long *I )
{
	int rank, offset, code, i, ks;
	/* - compute rank of permutation of x
     (x is a codevector in a permutation code C,
      the code C is organized according to the lexicographic order,
      the maximum xs of C is called the signed leader (playing the role of
      a generator for C),
      the rank #include <assert.h>of the permutation is the index of x in C)
    - get also the sign code of xs as a side product
      (the sign code is obtained by concatenating the sign bits of xs) */
	re8_compute_rank_of_permutation_and_sign_code( x, &rank, &code );
	/* compute cardinality offset in 2 steps:
     1. search for the sign code of xs in a pre-computed list of sign codes
        (the result ks is an identifier of the signed leader xs related to x,
	 the search is focused based on the identifier ka of the absolute
	 leader related to x)
     2. get the cardinality offset by table look-up
  */
	ks = -1; /* initialization to use assert() */
	for( i = Ia[*ka]; i < NB_LDSIGN; i++ ) {
		if( code == Ds[i] ) {
			ks = i;
			break;
		}
	}
	assert( ks >= 0 );
	offset = Is[ks];
	/* compute index of x: index = cardinality offset + rank of permutation */
	*I = offset + rank;
}
/*---------------------------------------------------------------------
  re8_compute_rank_and_sign_code(int *x, t, code)
  COMPUTE THE RANK OF THE PERMUTATION GIVEN x FROM THE SIGNED LEADER xs
  AND AS A SIDE-PRODUCT COMPUTE THE SIGN CODE OF xs
  (i) x:         point in RE8 (8-dimensional integer vector)
  (o) rank:      rank of the permutation x of xs (scalar integer)
  (o) sign_code: sign code of xs (8-bit word)
  note : the rank computation is based on Schalkwijk formula.
                 __
                \         (7-k)!
         rank =  |    ---------------   x (wk0 +...+ wk(d(k)-1))
                /__      ___
               k=0..7   |   |
                        |   |  (w_kî)!
                       i=0..q-1
   where k   =position in the vector
         q   =size of the alphabet of xs
         w_ki=number of occurences of the ith symbol from th kth position
	 d(k)=code translation a symbol of the alphabet into an index in
              the alphabet
  ---------------------------------------------------------------------
 */
void re8_compute_rank_of_permutation_and_sign_code( int *x, int *rank, int *sign_code )
{
	int xs[8], a[8], q, d[8], w[8], A, B, idx, tmp, abs_i, abs_j;
	int i, j, k;
	/* sort the elements of x to obtain the signed leader xs of x */
	for( i = 0; i < 8; i++ ) {
		xs[i] = x[i];
	}
	for( k = 0; k < 7; k++ ) {
		j = k;
		for( i = k + 1; i < 8; i++ ) {
			abs_j = abs( xs[j] );
			abs_i = abs( xs[i] );
			if( abs_i >= abs_j ) {
				if( abs_i > xs[j] ) {
					j = i;
				}
			}
		}
		if( j > k ) {
			tmp = xs[k];
			xs[k] = xs[j];
			xs[j] = tmp;
		}
	}
	/* compute the sign code of xs (this is a side-product of the rank
     calculation) */
	*sign_code = 0;
	for( i = 0; i < 8; i++ ) {
		if( xs[i] < 0 ) {
			*sign_code += tab_pow2[i]; /* *sign_code += 1<<(7-i); */
		}
	}
	/* compute the alphabet a=[a[0] ... a[q-1]] of x (q elements)
     such that a[0] != ... != a[q-1] and compute q */
	a[0] = xs[0];
	q = 1;
	for( i = 1; i < 8; i++ ) {
		if( xs[i] != xs[i - 1] ) {
			a[q] = xs[i];
			q++;
		}
	}
	/* translate x into d (where 0<=d[i]<q and d[i]=j if x[i]=a[j])
     based on the alphabet a */
	for( i = 0; i < 8; i++ ) {
		for( j = 0; j < q; j++ ) {
			if( x[i] == a[j] ) {
				d[i] = j;
				break; /* end loop over j */
			}
		}
	}
	/* compute rank of permutation based on Schalkwijk's formula
     the rank is given by rank=sum_{k=0..7} (A_k * fac_k/B_k) */
	*rank = 0;
	for( j = 0; j < q; j++ ) {
		w[j] = 0;
	}
	B = 1;
	for( i = 7; i >= 0; i-- ) {
		idx = d[i];
		w[idx]++;
		B *= w[idx];
		A = 0;
		for( j = 0; j < idx; j++ ) {
			A += w[j];
		}
		if( A > 0 ) {
			/* A_k * fac_k/B_k */
			*rank += A * tab_factorial[i] / B;
		}
	}
}
