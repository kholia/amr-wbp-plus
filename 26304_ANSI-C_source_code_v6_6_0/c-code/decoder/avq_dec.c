/*---------------------------------------------------------------------------*
 *         SPLIT ALGEBRAIC VECTOR QUANTIZER BASED ON RE8 LATTICE             *
 *---------------------------------------------------------------------------*
 * NOTE: a mitsmatch can occurs in some subvectors between the encoder       *
 *       and decoder, because the encoder use a bit-rate estimator to set    *
 *       the TCX global gain - this estimator is many times faster than the  *
 *       call of RE8_idx() for bits calculation.                             *
 *---------------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NQ_MAX 36
#define FAC_LOG2 3.321928095
#define NSV_MAX 256 /* number of sub-vector max in QVAE, 256*8=2048 */
/* local functions */
static void read_nq( int *nq, int *n_bits, int *parm, int *pos );
static void read_all_nq( int n_pack, int NB_BITS, int Nsv, int *pos_n, int *nq, int *last, int **parm_ptr );
static void chk_ovf( int n_bits, int n, int *n1, int *n2 );
static void read_I( int n, int *pos_i, long *I, int *parm );
static void read_k( int n, int *pos_i, int *k, int *parm, int flag );
static void init_pos_i_ovf( int n_pack, int *nq, int *pos_n, int last, int *pos_i_ovf );
void        read_ovf( int n_pack, int *parm_ovf, int n, int *pos_i_ovf, int *pos_n, int **parm_ptr );
static void read_all_i( int n_pack, int *nq, int *pos_n, int last, long *I, int *kv, int **parm_ptr );
static void read_track( int *pos_n, int *pos_i, int *nq, long *I, int *kv, int *parm );
/*
  AVQ_demuxdec_multi(parm, xriq, parm, NB_BITS, Nsv, n_pack, bfi)
  DEMULTIPLEX AND DECODE SUBVECTORS FROM SEVERAL PACKETS
  -> n_pack  : number of packets
  -> parm    : bitstream = table of 4-bit words [0..(NB_BITS-1)/4]
               the bitstream is divided into n_pack blocks
             parm[0]..parm[NB_BITS-1]
  <- xriq    : rounded subvectors [0..8*Nsv-1]
               followed by rounded bit allocations [8*Nsv..8*Nsv+Nsv-1]
  -> NB_BITS : number of bits allocated for split multi-rate RE8 VQ
  -> Nsv     : number of subvectors
  -> bfi     : bad frame indicator for each packet
  notes (see AVQ_encmux_multi)
*/
void AVQ_demuxdec( int n_pack, int *param, int *n_bits, float *xriq, int Nsv, int *bfi )
{
	int  p, pos_i[N_PACK_MAX], pos_n[N_PACK_MAX], l, i;
	int  any_loss;
	int  nq[NSV_MAX], *kv, last, c[8];
	long I[NSV_MAX];
	int *parm_ptr[4];
	int  NB_BITS;
	kv = (int *)xriq; /* reuse vector to save memory */
	NB_BITS = 0;
	for( i = 0; i < n_pack; i++ ) {
		parm_ptr[i] = param;
		NB_BITS += n_bits[i];
		param += ( ( n_bits[i] + 3 ) / 4 );
	}
	/* initialize pointers */
	for( p = 0; p < n_pack; p++ ) {
		pos_i[p] = 0;
		pos_n[p] = n_bits[p] - 1;
	}
	/* verify if any packet is lost */
	any_loss = 0;
	for( p = 0; p < n_pack; p++ ) {
		any_loss |= bfi[p];
	}
	/* demultiplex and decode */
	if( any_loss == 0 ) {
		/* read and decode encoded nq, set overflow pointers and compute
       number of unused bits in each packet and index of last described
       subvector */
		read_all_nq( n_pack, NB_BITS, Nsv, pos_n, nq, &last, parm_ptr );
		/* read indices (i) */
		read_all_i( n_pack, nq, pos_n, last, I, kv, parm_ptr );
		/* decode all */
		for( l = 0; l < Nsv; l++ ) {
			/* multi-rate RE8 decoder */
			RE8_dec( nq[l], I[l], &kv[8 * l], c );
			/* write decoded RE8 vector to decoded subvector #l */
			for( i = 0; i < 8; i++ ) {
				xriq[l * 8 + i] = (float)c[i];
			}
		}
	}
	else {
		/* do not decode bits in overflow if any packet is lost */
		for( p = 0; p < n_pack; p++ ) {
			if( bfi[p] == 0 ) {
				/* read and decode encoded codebook numbers and indices jointly
           stop when oveflow or after last subvector */
				for( l = p; l < Nsv; l += n_pack ) {
					/* read and decode codebook number (nq), read index (i)
             note : we cannot use read_parm because of pos_i/pos_n
                    which may not be multiple of 4 */
					read_track( pos_n + p, pos_i + p, &nq[l], &I[l], &kv[8 * l], parm_ptr[p] );
					/* multi-rate RE8 decoder */
					RE8_dec( nq[l], I[l], &kv[8 * l], c );
					/* write decoded RE8 vector to subvector l */
					for( i = 0; i < 8; i++ ) {
						xriq[l * 8 + i] = (float)c[i];
					}
				}
			}
			else {
				for( l = p; l < Nsv; l += n_pack ) {
					for( i = 0; i < 8; i++ ) {
						xriq[8 * l + i] = 0.0f;
					}
				}
			}
		}
	}
	return;
}
/* read a single codebook number */
static void read_nq( int *nq, int *n_bits, int *parm, int *pos )
{
	*nq = 0;
	if( *n_bits >= 9 ) {
		if( ( ( parm[*pos / 4] >> ( *pos % 4 ) ) & 0x01 ) == 1 ) {
			*nq = 2;
			*n_bits -= 9;
			*pos = *pos - 1;
			/* CHECK: add a test for nq == NQ_MAX */
			while( ( ( ( parm[*pos / 4] >> ( *pos % 4 ) ) & 0x01 ) == 1 ) && ( *n_bits >= 5 ) && ( *nq < NQ_MAX ) ) {
				*nq = *nq + 1;
				*n_bits -= 5;
				*pos = *pos - 1;
			}
		}
	}
	if( *n_bits > 0 ) {
		*n_bits = *n_bits - 1;
		*pos = *pos - 1;
	}
}
/*
  read_all_nq(pos_n, n_pack, nq, pos_ovf, n_bits_left, last, parm)
  DE-MULTIPLEX AND DECODE AND ALL CODEBOOKS NUMBERS IN nq[] TRACK-BY-TRACK
   ->  nq         : table of codebook numbers [0..Nsv-1]
   ->  pos_n      : table of pointers to write nq in packets [0..n_pack-1]
   ->  NB_BITS    : number of bits allocated for split multi-rate RE8 VQ
   ->  n_pack     : number of packets
   ->  last       : index of last described subvector
  <-  pos_ovf     : pointers for overflow [0..n_pack-1]
  <-  n_bits_left : number of unused bits in packets [0..n_pack-1]
  <-> parm        : bistream
*/
static void read_all_nq( int n_pack, int NB_BITS, int Nsv, int *pos_n, int *nq, int *last, int **parm_ptr )
{
	int n_bits, p, l;
	n_bits = NB_BITS;
	*last = -1;
	for( l = 0; l < Nsv; l++ ) {
		p = l % n_pack;
		read_nq( &nq[l], &n_bits, parm_ptr[p], &pos_n[p] );
		if( nq[l] > 0 ) {
			*last = l;
		}
	}
	return;
}
/* check if n groups of 4 bits (i.e. 4n bits) fit in n_bits bits */
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
/* read 4n bits for base codebook index (I) */
static void read_I( int n, int *pos_i, long *I, int *parm )
{
	int pos;
	/* base codebook index */
	*pos_i = *pos_i + 4 * n;
	pos = *pos_i / 4 - 1;
	while( n-- > 0 ) {
		*I = *I << 4;
		*I = *I + ( parm[pos--] & 0x0F );
	}
}
/* read 4n bits for Voronoi index (k[]) */
static void read_k( int n, int *pos_i, int *k, int *parm, int flag )
{
	int pos, i, ival, delta, *kv;
	*pos_i = *pos_i + 4 * n;
	pos = *pos_i / 4 - 1;
	delta = 4 * flag;
	while( n-- > 0 ) {
		kv = k + delta;
		ival = ( parm[pos--] & 0x0F ); /* optional mask */
		for( i = 3; i >= 0; i-- ) {
			kv[i] <<= 1;
			kv[i] += ( ival & 0x01 );
			ival >>= 1;
		}
		delta = ( delta + 4 ) % 8; /* circular shift */
	}
}
/* split a codebook number (nq) into a number of bits for the base
   codebook index (4 x ni) and  the Voronoi index (4 x nk) */
static void split_n( int nq, int *ni, int *nk )
{
	int tmp;
	*ni = nq;
	*nk = 0;
	if( *ni > 4 ) {
		tmp = ( *ni - 4 + 1 ) / 2;
		*nk = tmp * 2;
		*ni -= *nk;
	}
}
/* find in each packet the positions where overflow occurs */
static void init_pos_i_ovf( int n_pack, int *nq, int *pos_n, int last, int *pos_i_ovf )
{
	int p, pos, n_bits, l, n1, n2;
	/* find in each packet the positions where overflow occurs & count the number
     of bits to put in the extra packet */
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
/* read bits in overflow */
void read_ovf( int n_pack, int *parm_ovf, int n, int *pos_i_ovf, int *pos_n, int **parm_ptr )
{
	int p, n_bits, pos, *parm, pos_ovf, moved_bit;
	/* initialize position in overflow packet (parm_ovf[]) */
	pos_ovf = 0;
	/* read 4-bit chunks */
	for( p = 0; p < n_pack; p++ ) {
		n_bits = pos_n[p] - pos_i_ovf[p] + 1;
		if( ( n_bits >= 4 ) && ( n > 0 ) ) {
			pos = pos_i_ovf[p] / 4;
			parm = parm_ptr[p];
			do {
				parm_ovf[pos_ovf++] = parm[pos++];
				n_bits -= 4;
				n -= 4;
			} while( ( n_bits >= 4 ) && ( n > 0 ) );
			pos_i_ovf[p] = pos * 4;
		}
	}
	pos_ovf *= 4;
	/* read bit-by-bit */
	for( p = 0; p < n_pack; p++ ) {
		n_bits = pos_n[p] - pos_i_ovf[p] + 1;
		if( ( n_bits > 0 ) && ( n > 0 ) ) {
			pos = pos_i_ovf[p];
			parm = parm_ptr[p];
			do {
				moved_bit = ( parm[pos / 4] >> ( pos % 4 ) ) & 0x01;
				parm_ovf[pos_ovf / 4] += moved_bit << ( pos_ovf % 4 );
				n_bits--;
				n--;
				pos++;
				pos_ovf++;
			} while( ( n_bits > 0 ) && ( n > 0 ) );
			pos_i_ovf[p] = pos;
		}
	}
}
/*
  read_all_i(nq, p, pos_n, n_pack, last, pos_ovf, n_bits_left, I, kv, parm)
  DEMULTIPLEX AN INDEX (I,kv) OF 4*nq BITS FROM BITSTREAM
  (THE INDEX CORRESPONDS TO A SUBVECTOR IN TRACK #p)
   -> n_pack      : number for packets
   -> nq          : codebook numbers
   -> pos_n       : position in bitstream after reading nq[]
   -> last        : index of last described subvector
  <-  I           : base codebook index
  <-  kv          : Voronoi index
   -> parm_ptr    : multiple bistreams
*/
static void read_all_i( int n_pack, int *nq, int *pos_n, int last, long *I, int *kv, int **parm_ptr )
{
	int  pos_ovf, p, l, n_bits, ni, nk, n1, i, n2, pos, *parm, pos_i_ovf[N_PACK_MAX];
	long index;
	int  parm_ovf[NQ_MAX];
	/* initialize overflow packet */
	for( i = 0; i < NQ_MAX; i++ ) {
		parm_ovf[i] = 0;
	}
	/* find positions in each packet from which the bits in overflow can be read */
	init_pos_i_ovf( n_pack, nq, pos_n, last, pos_i_ovf );
	/* read indices */
	for( p = 0; p < n_pack; p++ ) {
		pos = 0; /* pos_i[p]/4 */
		parm = parm_ptr[p];
		for( l = p; l <= last; l += n_pack ) {
			if( nq[l] > 0 ) {
				/* compute number of bits left for indices in packet #p */
				n_bits = pos_n[p] - pos + 1;
				split_n( nq[l], &ni, &nk );
				/* read I in packet #p and in overflow */
				chk_ovf( n_bits, ni, &n1, &n2 );
				index = 0;
				if( n2 > 0 ) {
					pos_ovf = 0;
					read_ovf( n_pack, parm_ovf, 4 * n2, pos_i_ovf, pos_n, parm_ptr );
					read_I( n2, &pos_ovf, &index, parm_ovf );
					for( i = 0; i < n2; i++ ) {
						parm_ovf[i] = 0;
					}
				}
				read_I( n1, &pos, &index, parm );
				n_bits -= 4 * n1;
				I[l] = index;
				/* read Voronoi index  */
				if( nk > 0 ) {
					for( i = 0; i < 8; i++ ) {
						kv[8 * l + i] = 0;
					}
					chk_ovf( n_bits, nk, &n1, &n2 );
					if( n2 > 0 ) {
						pos_ovf = 0;
						read_ovf( n_pack, parm_ovf, 4 * n2, pos_i_ovf, pos_n, parm_ptr );
						read_k( n2, &pos_ovf, &kv[8 * l], parm_ovf, 1 );
						for( i = 0; i < n2; i++ ) {
							parm_ovf[i] = 0;
						}
					}
					read_k( n1, &pos, &kv[8 * l], parm, ( n2 + 1 ) % 2 );
				}
			}
		}
	}
}
/*
  read_track(pos_n, pos_i, nq, I, kv, parm)
  DEMULTIPLEX nq AND INDEX (I,kv) OF 4*nq BITS FROM BITSTREAM
  (THE INDEX CORRESPONDS TO A SUBVECTOR IN TRACK #p)
  <-> pos_n       : pointer to read nq
  <-> pos_i       : pointer to read i
  <-  nq          : codebook number (scalar)
  <-  I           : base codebook index
  <-  kv          : Voronoi index
   -> parm        : bistream
*/
static void read_track( int *pos_n, int *pos_i, int *nq, long *I, int *kv, int *parm )
{
	int n_bits, i, ni, nk;
	/* compute number of bits left for indices in packet #p */
	n_bits = *pos_n - *pos_i + 1;
	/* read nq */
	read_nq( nq, &n_bits, parm, pos_n );
	/* read i and kv */
	if( *nq > 0 ) {
		*I = 0;
		for( i = 0; i < 8; i++ ) {
			kv[i] = 0;
		}
		split_n( *nq, &ni, &nk );
		read_I( ni, pos_i, I, parm );
		read_k( nk, pos_i, kv, parm, 1 );
	}
}
