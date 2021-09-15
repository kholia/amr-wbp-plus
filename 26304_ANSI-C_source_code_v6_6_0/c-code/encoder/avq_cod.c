/*---------------------------------------------------------------------------*
 *         SPLIT ALGEBRAIC VECTOR QUANTIZER BASED ON RE8 LATTICE             *
 *---------------------------------------------------------------------------*
 * NOTE: a mitsmatch can occurs in some subvectors between the encoder       *
 *       and decoder, because the encoder use a bit-rate estimator to set    *
 *       the TCX global gain - this estimator is many times faster than the  *
 *       call of RE8_idx() for bits calculation.                             *
 *---------------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NQ_MAX 36
#define FAC_LOG2 3.321928095f
#define NSV_MAX 256 /* number of sub-vector max in QVAE, 256*8=2048 */
/* local function */
static int  calc_bits( int nq );
static void sort( int *ebits, int n, int *idx );
static void split_idx_noovf( int *xriq, int NB_BITS, int Nsv, int *nq, long *I, int *kv );
static void writ_all_nq( int n_pack, int *nq, int *pos_n, int NB_BITS, int Nsv, int *last, int **parm_ptr );
static void chk_ovf( int n_bits, int n, int *n1, int *n2 );
static void writ_I( int n, int *pos_i, long *I, int *parm );
static void writ_k( int n, int *pos_i, int *k, int *parm, int flag );
static void init_pos_i_ovf( int n_pack, int *nq, int *pos_n, int last, int *pos_i_ovf );
static void writ_ovf( int n_pack, int *parm_ovf, int n, int *pos_i_ovf, int *pos_n, int **parm_ptr );
static void writ_all_i( int n_pack, int *nq, int *pos_n, int last, long *I, int *kv, int **parm_ptr );
float       AVQ_cod(      /* output: comfort noise gain factor      */
    float *xri,     /* input:  vector to quantize             */
    int *  xriq,    /* output: quantized normalized vector (assuming the bit budget is enough) */
    int    NB_BITS, /* input:  number of allocated bits          */
    int    Nsv )       /* input:  number of subvectors (lg=Nsv*8) */
{
	int   i, l, n, iter, c[8];
	float gain, gain_inv, x1[8], ener, tmp, nbits, nbits_max, fac, offset;
	float ebits[NSV_MAX];
	/* find energy of each subvector in log domain (scaled for bits estimation) */
	for( l = 0; l < Nsv; l++ ) {
		for( i = 0; i < 8; i++ ) {
			x1[i] = xri[l * 8 + i];
		}
		ener = 2.0; /* to set ebits >= 0 */
		for( i = 0; i < 8; i++ ) {
			ener += x1[i] * x1[i];
		}
		/* estimated bit consumption when gain=1 */
		ebits[l] = 5.0f * FAC_LOG2 * (float)log10( ener * 0.5 );
	}
	/*---------------------------------------------------------------------*
   * subvector energy worst case:                                        *
   * - typically, it's a tone with maximum of amplitude (RMS=23170).     *
   * - fft length max = 1024 (N/2 is 512)                                *
   * log10(energy) = log10(23710*23710*1024*(N/2)) = 14.45               *
   * ebits --> 5.0*FAC_LOG2*14.45 = 240 bits                             *
   *---------------------------------------------------------------------*/
	/* estimate gain according to number of bits allowed */
	fac = 128.0; /* start at the middle (offset range = 0 to 255.75) */
	offset = 0.0;
	nbits_max = 0.95f * ( (float)( NB_BITS - Nsv ) );
	/* tree search with 10 iterations : offset with step of 0.25 bits (0.3 dB) */
	for( iter = 0; iter < 10; iter++ ) {
		offset += fac;
		/* calculate the required number of bits */
		nbits = 0.0;
		for( l = 0; l < Nsv; l++ ) {
			tmp = ebits[l] - offset;
			if( tmp < 0.0 ) {
				tmp = 0.0;
			}
			nbits += tmp;
		}
		/* decrease gain when no overflow occurs */
		if( nbits <= nbits_max ) {
			offset -= fac;
		}
		fac *= 0.5;
	}
	/* estimated gain (when offset=0, estimated gain=1) */
	gain = (float)pow( 10.0, offset / ( 2.0 * 5.0 * FAC_LOG2 ) );
	gain_inv = 1.0f / gain;
	/* quantize all subvector using estimated gain */
	for( l = 0; l < Nsv; l++ ) {
		for( i = 0; i < 8; i++ ) {
			x1[i] = xri[l * 8 + i] * gain_inv;
		}
		RE8_PPV( x1, c );
		for( i = 0; i < 8; i++ ) {
			xriq[l * 8 + i] = c[i];
		}
	}
	/* evaluate comfort noise level at freq. over 3200Hz (Nsv/2) */
	/* SV with ebits < 5 bits may not be quantized */
	nbits = 0;
	n = 1;
#ifdef ENHANCE_TCX_NOISE_FILL
	/* reduce noise level when SNR is high over 3200Hz (on music harmonic tone) */
	tmp = 0.0f;
	for( l = Nsv / 2; l < Nsv; l++ )
		if( ebits[l] > tmp )
			tmp = ebits[l];
	tmp = tmp - 10.0f;
	if( offset < tmp )
		offset = tmp;
#endif
	for( l = Nsv / 2; l < Nsv; l++ ) {
		tmp = ebits[l] - offset;
		if( tmp < 5.0 ) {
			nbits += tmp;
			n++;
		}
	}
	nbits /= (float)n;
	fac = (float)pow( 10.0, ( nbits - 5.0 ) / ( 2.0 * 5.0 * FAC_LOG2 ) );
	/* round bit allocations and save */
	for( i = 0; i < Nsv; i++ ) {
		xriq[( Nsv * 8 ) + i] = (int)floor( ebits[i] * 128.0 );
	}
	return ( fac );
}
/*
  AVQ_encmux(xriq, parm, n_bits, NB_BITS, Nsv, n_pack)
  ENCODE AND MULTIPLEX SUBVECTORS INTO SEVERAL PACKETS
   -> n_pack  : number of packets
   -> xriq    : rounded subvectors [0..8*Nsv-1]
                followed by rounded bit allocations [8*Nsv..8*Nsv+Nsv-1]
  <-> param   : multiplexed parameters
   -> n_bits  : size of each packet
   -> Nsv     : number of subvectors
  note:
  Nsv MUST be multiple of n_pack
  IMPORTANT:
  it is assumed that codebook numbers in track #p do not cause bit
  budget overflow in packet #p
  in practice this is ok if p<5 because the quantizer #n takes 5n bits
  and putting all bits in subvectors of track #p results in NB_BITS/5 bits
  for codebook numbers
*/
void AVQ_encmux( int n_pack, int *xriq, int *param, int *n_bits, int Nsv )
{
	int  last, i, p;
	int  kv[NSV_MAX * 8], nq[NSV_MAX];
	long I[NSV_MAX];
	int  pos_n[N_PACK_MAX];
	/* int pos_i[N_PACK_MAX]; */
	int *parm;
	int *parm_ptr[4];
	int  NB_BITS;
	NB_BITS = 0;
	for( i = 0; i < n_pack; i++ ) {
		parm_ptr[i] = param;
		NB_BITS += n_bits[i];
		param += ( ( n_bits[i] + 3 ) / 4 );
	}
	/* initialize pointers and packets */
	for( p = 0; p < n_pack; p++ ) {
		/* pos_i[p] = 0; */
		pos_n[p] = n_bits[p] - 1;
		/* initialize packet to zero, i.e. all subvectors are set to zero */
		parm = parm_ptr[p];
		for( i = 0; i <= ( ( n_bits[p] - 1 ) / 4 ); i++ ) {
			parm[i] = 0;
		}
	}
	/* encode subvectors and fix possible overflows in TOTAL bit budget:
     i.e. find (i,nq) for each subvector where
          i is a codevector index split in a base codebook index (I) and
                                                  a Voronoi index (kv)
          nq is a codebook index (nq=0,2,3,4,...) */
	split_idx_noovf( xriq, NB_BITS, Nsv, nq, I, kv );
	/* split multiplexing of codebook numbers (by interleaved tracks) */
	writ_all_nq( n_pack, nq, pos_n, NB_BITS, Nsv, &last, parm_ptr );
	/* write indices
     multiplexing is done track-by-track (from track #0 to track #n_pack-1) */
	writ_all_i( n_pack, nq, pos_n, last, I, kv, parm_ptr );
	return;
}
/*
  calc_bits(nq)
  COMPUTE (NUMBER OF BITS -1) TO DESCRIBE Q #nq
  -> nq: quantizer id (0,2,3,4...)
  <-   : bit allocation
*/
static int calc_bits( int nq )
{
	if( nq >= 2 ) {
		/* 4n bits + variable-length descriptor for allocation:
       descriptor -> nq
       0          -> 0
       10         -> 2
       110        -> 3
       => size of descriptor = 5n bits */
		return ( ( nq * 5 ) - 1 ); /* [5n-1] */
	}
	else {
		return 0; /* 1-1 [1 bit to describe the allocation] */
	}
}
/*
  sort(ebits, n, idx)
  SORT SUBVECTORS BY DECREASING BIT ALLOCATIONS
  -> ebits : estimated bit allocations (table of n *positive* integers)
  -> n     : number of subvectors
  <- idx   : indices
*/
static void sort( int *ebits, int n, int *idx )
{
	int t[NSV_MAX], i, j, ebits_max, pos;
	for( i = 0; i < n; i++ ) {
		t[i] = ebits[i];
	}
	for( i = 0; i < n; i++ ) {
		ebits_max = t[0];
		pos = 0;
		for( j = 1; j < n; j++ ) {
			if( t[j] > ebits_max ) {
				ebits_max = t[j];
				pos = j;
			}
		}
		idx[i] = pos;
		t[pos] = -1;
	}
	return;
}
/*
  split_idx_noovf(xriq,NB_BITS, Nsv, nq, I, kv, last)
  COMPUTE MULTI-RATE INDICES FOR ALL SUBVECTORS AND FORCE NO BIT BUDGET OVERFLOW
  PRIOR TO MULTIPLEXING
  -> xriq    : rounded subvectors [0..8*Nsv-1]
               followed by rounded bit allocations [8*Nsv..8*Nsv+Nsv-1]
  -> NB_BITS : number of bits allocated for split multi-rate RE8 VQ
  -> Nsv     : number of subvectors
  <- nq      : codebook numbers
  <- I       : indices for base quantizers (Q2,Q3,Q4)
  <- kv      : Voronoi indices
*/
static void split_idx_noovf( int *xriq, int NB_BITS, int Nsv, int *nq, long *I, int *kv )
{
	int k, l, n, n_bits, pos, pos_max;
	int sort_idx[NSV_MAX];
	/* sort subvectors by estimated bit allocations in decreasing order
     (l=idx[0] is such that (rounded) ebits[l] is maximum) */
	sort( &xriq[8 * Nsv], Nsv, sort_idx );
	/* compute multi-rate indices and avoid bit budget overflow  */
	pos_max = 0;
	n_bits = 0;
	for( l = 0; l < Nsv; l++ ) {
		/* find vector to quantize (criteria: nb of estimated bits) */
		pos = sort_idx[l];
		/* compute multi-rate index of rounded subvector (nq,I,kv[]) */
		RE8_cod( &xriq[pos * 8], &nq[pos], &I[pos], &kv[8 * pos] );
		if( nq[pos] > 0 ) {
			k = pos_max;
			if( pos > k ) {
				k = pos;
			}
			/* check for overflow and compute number of bits-1 (n) */
			n = calc_bits( nq[pos] );
			if( ( n_bits + n + k ) > NB_BITS ) { /* if budget overflow */
				nq[pos] = 0;                     /* force Q0 */
			}
			else {
				n_bits += n;
				pos_max = k; /* update index of last described subvector (last) */
			}
		}
	}
}
/*
  writ_all_nq(n_pack, nq, pos_n, NB_BITS, last, parm_ptr)
  ENCODE AND MULTIPLEX ALL CODEBOOKS NUMBERS IN nq[] TRACK-BY-TRACK
   ->  n_pack      : number of packets
   ->  nq          : table of codebook numbers [0..Nsv-1]
   ->  pos_n       : table of pointers to write nq in packets [0..n_pack-1]
   ->  NB_BITS     : total bit allocation
  <-  last        : index of last subvector for which an index is written
  <-> parm        : bistream
*/
static void writ_all_nq( int n_pack, int *nq, int *pos_n, int NB_BITS, int Nsv, int *last, int **parm_ptr )
{
	int p, pos, l, i, n_bits, *parm;
	n_bits = NB_BITS;
	*last = -1;
	/* write nq[l] for l=0...Nsv-1 in packet #p= mod(l,number_of_packets)*/
	for( l = 0; l < Nsv; l++ ) {
		p = l % n_pack; /* alternative : apply mask 0x01 / 0x03 */
		parm = parm_ptr[p];
		pos = pos_n[p];
		/* compute the minimal bit budget to write nq[l] and the related index */
		i = calc_bits( nq[l] );
		/* if the budget is exceeded, force nq[l] to 0
       else decrement the number of left bits */
		if( i > n_bits ) {
			nq[l] = 0;
		}
		else {
			n_bits -= i;
		}
		/* update "last" */
		if( nq[l] >= 2 ) {
			*last = l;
		}
		/* write the unary code (except stop bit) for nq[l] in packet #p */
		i = nq[l] - 1;
		while( i-- > 0 ) {
			parm[pos / 4] += 1 << ( pos % 4 );
			pos--;
		}
		/* if bit budget is not empty, write stop bit of unary code */
		if( n_bits > 0 ) {
			pos--;
			n_bits--;
		}
		pos_n[p] = pos;
	}
	return;
}
/* check if n groups of 4 bits fit in n_bits bits */
static void chk_ovf( int n_bits, int n, int *n1, int *n2 )
{
	if( 4 * n <= n_bits ) {
		*n1 = n;
		*n2 = 0;
	}
	else {
		*n1 = n_bits / 4; /* >> 2*/
		*n2 = n - *n1;
	}
}
/* write n groups of 4-bit for base codebook index (I) */
static void writ_I( int n, int *pos_i, long *I, int *parm )
{
	int pos;
	/* base codebook index */
	pos = *pos_i / 4;
	while( n-- > 0 ) {
		parm[pos++] = ( *I & 0x0F );
		*I = *I >> 4;
	}
	*pos_i = pos * 4;
}
/* write n groups of 4-bit for Voronoi index (k[]) */
static void writ_k( int n, int *pos_i, int *k, int *parm, int flag )
{
	int i, ival, delta, *kv, pos;
	delta = 4 * flag;
	pos = *pos_i / 4;
	while( n-- > 0 ) {
		kv = k + delta;
		ival = 0;
		for( i = 0; i < 4; i++ ) {
			ival <<= 1;
			ival += ( kv[i] & 0x01 );
			kv[i] >>= 1;
		}
		parm[pos++] = ival;
		delta = ( delta + 4 ) % 8; /* circular shift */
	}
	*pos_i = pos * 4;
}
/* find in each packet the positions where overflow occurs */
static void init_pos_i_ovf( int n_pack, int *nq, int *pos_n, int last, int *pos_i_ovf )
{
	int p, pos, n_bits, l, n1, n2;
	for( p = 0; p < n_pack; p++ ) {
		pos = 0;
		n_bits = pos_n[p] + 1; /* pos_n[p] - 0 +1 */
		for( l = p; l <= last; l += n_pack ) {
			if( nq[l] > 0 ) {
				chk_ovf( n_bits, nq[l], &n1, &n2 );
				n_bits -= 4 * n1;
				pos += n1;
			}
		}
		pos_i_ovf[p] = pos * 4;
	}
}
/* write bits in overflow */
static void writ_ovf( int n_pack, int *parm_ovf, int n, int *pos_i_ovf, int *pos_n, int **parm_ptr )
{
	int pos_ovf, p, n_bits, pos, *parm, moved_bit;
	/* initialize position in overflow packet (parm_ovf[]) */
	pos_ovf = 0;
	/* move bits from overflow packet (parm_ovf[]) to packets (parm_ptr[][])
     [write 4-bit by 4-bit] */
	for( p = 0; p < n_pack; p++ ) {
		/* compute number of bits left in packet #p */
		n_bits = pos_n[p] - pos_i_ovf[p] + 1;
		if( ( n_bits >= 4 ) && ( n > 0 ) ) {
			pos = pos_i_ovf[p] / 4;
			parm = parm_ptr[p];
			do {
				parm[pos++] = parm_ovf[pos_ovf++]; /* move 4 bits */
				n_bits -= 4;
				n -= 4;
			} while( ( n_bits >= 4 ) && ( n > 0 ) );
			pos_i_ovf[p] = pos * 4;
		}
	}
	pos_ovf *= 4;
	/* move bits remaining in overflow packet
     [write bit-by-bit (3 bits at maximum per packet)] */
	for( p = 0; p < n_pack; p++ ) {
		/* compute number of bits left in packet #p */
		n_bits = pos_n[p] - pos_i_ovf[p] + 1;
		if( ( n_bits > 0 ) && ( n > 0 ) ) {
			pos = pos_i_ovf[p];
			parm = parm_ptr[p];
			do {
				/* write a single bit */
				moved_bit = ( parm_ovf[pos_ovf / 4] >> ( pos_ovf % 4 ) ) & 0x01;
				parm[pos / 4] += moved_bit << ( pos % 4 );
				pos++;
				pos_ovf++;
				n_bits--;
				n--;
			} while( ( n_bits > 0 ) && ( n > 0 ) );
			pos_i_ovf[p] = pos;
		}
	}
}
/*
  writ_all_i(nq, p, pos_n, n_pack, last, pos_ovf, n_bits_left, parm_ptr)
  MULTIPLEX AN INDEX (I,kv) OF 4*nq BITS INTO BITSTREAM
  (THE INDEX CORRESPONDS TO A SUBVECTOR IN TRACK #p)
  ->  nq          : codebook number (scalar)
  ->  p           : index of track
  ->  pos_i       : pointer to write i in packet [0..n_pack-1]
  ->  n_pack      : number of packets
  <-> pos_ovf     : pointers for overflow [0..n_pack-1]
  <-> n_bits_left : number of unused bits in packets [0..n_pack-1]
  ->  I           : base codebook index
  ->  kv          : Voronoi index
  <-> parm        : bistream [initialized to zero when writ_all_i is called]
  important note:
  - if the index fits completely in packet #p, multiplexing is done as
    in writ_parm from top to bottom:
    the 1st bit of I is written at the top
    the last bit of kv at the bottom
*/
static void writ_all_i( int n_pack, int *nq, int *pos_n, int last, long *I, int *kv, int **parm_ptr )
{
	int  pos_ovf, p, l, n_bits, ni, nk, n1, i, n2, pos, k[8], *parm;
	long index;
	int  parm_ovf[NQ_MAX];
	int  pos_i_ovf[N_PACK_MAX];
	/* initialize overflow packet */
	for( i = 0; i < NQ_MAX; i++ ) {
		parm_ovf[i] = 0;
	}
	/* find positions in each packet from which the bits in overflow can be
     written */
	init_pos_i_ovf( n_pack, nq, pos_n, last, pos_i_ovf );
	/* write indices */
	for( p = 0; p < n_pack; p++ ) {
		pos = 0;
		parm = parm_ptr[p];
		for( l = p; l <= last; l += n_pack ) {
			if( nq[l] > 0 ) {
				/* compute number of bits left in packet #p */
				n_bits = pos_n[p] - pos + 1;
				/* compute number of 4-bit groups for base codebook index (ni)
           and Voronoi index (nk) */
				ni = nq[l];
				nk = 0;
				if( ni > 4 ) {
					nk = ( ni - 4 + 1 ) / 2; /* nkv*2 = number of 4-bit groups */
					nk *= 2;
					ni -= nk;
				}
				/* write base codebook index (in packet #p / overflow packet) */
				index = I[l];
				chk_ovf( n_bits, ni, &n1, &n2 );
				writ_I( n1, &pos, &index, parm );
				n_bits -= 4 * n1;
				if( n2 > 0 ) {
					/* write 4-bit groups in overflow packet */
					pos_ovf = 0;
					writ_I( n2, &pos_ovf, &index, parm_ovf );
					/* distribute bits from overflow packet to packets
             #0 to n_pack-1 */
					writ_ovf( n_pack, parm_ovf, 4 * n2, pos_i_ovf, pos_n, parm_ptr );
					for( i = 0; i < n2; i++ ) {
						parm_ovf[i] = 0;
					}
				}
				if( nk > 0 ) {
					/* write Voronoi index (in packet #p / overflow packet) */
					for( i = 0; i < 8; i++ ) {
						k[i] = kv[8 * l + i];
					}
					chk_ovf( n_bits, nk, &n1, &n2 );
					writ_k( n1, &pos, k, parm, 0 );
					if( n2 > 0 ) {
						pos_ovf = 0;
						writ_k( n2, &pos_ovf, k, parm_ovf, n1 % 2 );
						writ_ovf( n_pack, parm_ovf, 4 * n2, pos_i_ovf, pos_n, parm_ptr );
						for( i = 0; i < n2; i++ ) {
							parm_ovf[i] = 0;
						}
					}
				}
			}
		}
	}
}
