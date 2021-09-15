#include "../include/amr_plus.h"
#include <assert.h>
#include <math.h>
extern void re8_k2y( int *k, int m, int *y );
static void re8_decode_base_index( int *n, long *I, int *y );
static void re8_decode_rank_of_permutation( int t, int *xs, int *x );
/*--------------------------------------------------------------------------
  RE8_dec(n, I, k, y)
  MULTI-RATE INDEXING OF A POINT y in THE LATTICE RE8 (INDEX DECODING)
  (i) n: codebook number (*n is an integer defined in {0,2,3,4,..,n_max})
  (i) I: index of c (pointer to unsigned 16-bit word)
  (i) k: index of v (8-dimensional vector of binary indices) = Voronoi index
  (o) y: point in RE8 (8-dimensional integer vector)
  note: the index I is defined as a 32-bit word, but only
  16 bits are required (long can be replaced by unsigned integer)
  --------------------------------------------------------------------------
 */
void RE8_dec( int n, long I, int k[], int y[] )
{
	int i, m, v[8];
	/* decode the sub-indices I and kv[] according to the codebook number n:
     if n=0,2,3,4, decode I (no Voronoi extension)
     if n>4, Voronoi extension is used, decode I and kv[] */
	if( n <= 4 ) {
		re8_decode_base_index( &n, &I, y );
	}
	else {
		/* compute the Voronoi modulo m = 2^r where r is extension order */
		m = 1;
		while( n > 4 ) {
			m *= 2;
			n -= 2;
		}
		/* decode base codebook index I into c (c is an element of Q3 or Q4)
       [here c is stored in y to save memory] */
		re8_decode_base_index( &n, &I, y );
		/* decode Voronoi index k[] into v */
		re8_k2y( k, m, v );
		/* reconstruct y as y = m c + v (with m=2^r, r integer >=1) */
		for( i = 0; i < 8; i++ ) {
			y[i] = m * y[i] + v[i];
		}
	}
	return;
}
/*--------------------------------------------------------------------------
  re8_decode_base_index(n, I, y)
  DECODING OF AN INDEX IN Qn (n=0,2,3 or 4)
  (i) n: codebook number (*n is an integer defined in {0,2,3,4})
  (i) I: index of c (pointer to unsigned 16-bit word)
  (o) y: point in RE8 (8-dimensional integer vector)
  note: the index I is defined as a 32-bit word, but only
  16 bits are required (long can be replaced by unsigned integer)
  --------------------------------------------------------------------------
 */
static void re8_decode_base_index( int *n, long *I, int *y )
{
	int  i, im, t, sign_code, ka, ks, rank, leader[8];
	long offset;
	if( ( *n == 4 ) && ( *I > 65519L ) ) {
		*I = 0;
	}
	if( *n < 2 ) {
		for( i = 0; i < 8; i++ ) {
			y[i] = 0;
		}
	}
	else {
		/* search for the identifier ka of the absolute leader (table-lookup)
       Q2 is a subset of Q3 - the two cases are considered in the same branch
     */
		/* switch <=> if (n==4)... else ... end */
		switch( *n ) {
		case 2:
		case 3:
			for( i = 1; i < NB_LDQ3; i++ ) {
				if( *I < (long)I3[i] )
					break;
			}
			ka = A3[i - 1];
			break;
		case 4:
			for( i = 1; i < NB_LDQ4; i++ ) {
				if( *I < (long)I4[i] )
					break;
			}
			ka = A4[i - 1];
			break;
		}
		/* reconstruct the absolute leader */
		/* Da[8*ka+i] -> pointer to Da[ka<<3] =>*/
		for( i = 0; i < 8; i++ ) {
			leader[i] = Da[ka][i];
		}
		/* search for the identifier ks of the signed leader (table look-up)
       (this search is focused based on the identifier ka of the absolute
        leader)*/
		t = Ia[ka];
		im = Ns[ka];
		for( i = im - 1; i >= 0; i-- ) {
			if( *I >= (long)Is[t + i] ) {
				ks = i;
				break;
			}
		}
		/* reconstruct the signed leader from its sign code */
		sign_code = Ds[t + ks];
		for( i = 0; i < 8; i++ ) {
			if( sign_code >= tab_pow2[i] ) {
				sign_code -= tab_pow2[i];
				leader[i] = -leader[i];
			}
		}
		/* compute the cardinality offset */
		offset = Is[t + ks];
		/* compute and decode the rank of the permutation */
		rank = *I - offset;
		re8_decode_rank_of_permutation( rank, leader, y );
	}
	return;
}
/*--------------------------------------------------------------------------
  re8_decode_rank_of_permutation(rank, xs, x)
  DECODING OF THE RANK OF THE PERMUTATION OF xs
  (i) rank: index (rank) of a permutation
  (i) xs:   signed leader in RE8 (8-dimensional integer vector)
  (o) x:    point in RE8 (8-dimensional integer vector)
  --------------------------------------------------------------------------
 */
static void re8_decode_rank_of_permutation( int rank, int *xs, int *x )
{
	int        a[8], q, w[8], B, A, fac, *ptr_w, *ptr_a;
	const int *ptr_factorial;
	int        i, j;
	long       target;
	/* --- pre-processing based on the signed leader xs ---
     - compute the alphabet a=[a[0] ... a[q-1]] of x (q elements)
       such that a[0]!=...!=a[q-1]
       it is assumed that xs is sorted in the form of a signed leader
       which can be summarized in 2 requirements:
          a) |xs[0]| >= |xs[1]| >= |xs[2]| >= ... >= |xs[7]|
          b) if |xs[i]|=|xs[i-1]|, xs[i]>=xs[i+1]
       where |.| indicates the absolute value operator
     - compute q (the number of symbols in the alphabet)
     - compute w[0..q-1] where w[j] counts the number of occurences of
       the symbol a[j] in xs
     - compute B = prod_j=0..q-1 (w[j]!) where .! is the factorial */
	/* xs[i], xs[i-1] and ptr_w/a*/
	ptr_w = w;
	ptr_a = a;
	*ptr_w = 1;
	*ptr_a = xs[0];
	q = 1;
	B = 1;
	for( i = 1; i < 8; i++ ) {
		if( xs[i] != xs[i - 1] ) {
			ptr_w++;
			ptr_a++;
			*ptr_w = 0;
			*ptr_a = xs[i];
			q++;
		}
		( *ptr_w )++;
		B *= *ptr_w;
	}
	/* --- actual rank decoding ---
     the rank of x (where x is a permutation of xs) is based on
     Schalkwijk's formula
     it is given by rank=sum_{k=0..7} (A_k * fac_k/B_k)
     the decoding of this rank is sequential and reconstructs x[0..7]
     element by element from x[0] to x[7]
     [the tricky part is the inference of A_k for each k...]
   */
	/* decode x element by element */
	ptr_factorial = tab_factorial;
	for( i = 0; i < 8; i++ ) {
		/* infere A (A_k): search j such that x[i] = a[j]
	 A = sum_{i=0...j-1} w[j] with 0<=j<q and sum_{i=0..-1} . = 0
	 j can be found by accumulating w[j] until the sum is superior to A
         [note: no division in search for j, but 32-bit arithmetic required] */
		fac = *ptr_factorial;
		target = -(long)rank * (long)B;
		A = 0;
		j = 0;
		do {
			target += w[j] * fac;
			A += w[j] * fac;
			j++;
		} while( target <= 0 );
		j--;
		assert( j < q );
		A -= w[j] * fac;
		x[i] = a[j];
		/* update rank, denominator B (B_k) and counter w[j] */
		if( A > 0 ) {
			rank -= A / B;
		}
		if( w[j] > 1 ) {
			B = B / w[j];
		}
		w[j]--;
		ptr_factorial++;
	}
}
