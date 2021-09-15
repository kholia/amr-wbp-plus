#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
static int unpack4bits( int nbits, int *prm, short *ptr );
/*-----------------------------------------------------------------*
 * Funtion  enc_prm()                                              *
 * ~~~~~~~~~~~~~~~~~~                                              *
 * encode AMR-WB+ parameters according to selected mode            *
 *-----------------------------------------------------------------*/
void enc_prm(
    int   mod[],      /* (i) : frame mode (mode[4], 4 division) */
    int   codec_mode, /* (i) : AMR-WB+ mode (see cnst.h)        */
    int   param[],    /* (i) : parameters                       */
    short serial[],   /* (o) : serial bits stream               */
    int   nbits_pack  /* (i) : number of bits per packet of 20ms*/
)
{
	int    j, k, n, sfr, mode, n_pack, nbits, *prm, parity, bit;
	int    nbits_AVQ[N_PACK_MAX], prm_AVQ[( NBITS_MAX / 4 ) + N_PACK_MAX];
	short *ptr;
	/* remove bits for mode (2 bits per 20ms packet) */
	nbits = ( NBITS_CORE[codec_mode] / 4 ) - 2;
	k = 0;
	while( k < NB_DIV ) {
		mode = mod[k];
		/* set pointer to parameters */
		prm = param + ( k * NPRM_DIV );
		if( ( mode == 0 ) || ( mode == 1 ) ) {
			ptr = serial + ( k * nbits_pack ) + 2; /* +2 because of mode */
			/* encode LPC parameters (46 bits) */
			int2bin( prm[0], 8, ptr );
			ptr += 8;
			int2bin( prm[1], 8, ptr );
			ptr += 8;
			int2bin( prm[2], 6, ptr );
			ptr += 6;
			int2bin( prm[3], 7, ptr );
			ptr += 7;
			int2bin( prm[4], 7, ptr );
			ptr += 7;
			int2bin( prm[5], 5, ptr );
			ptr += 5;
			int2bin( prm[6], 5, ptr );
			ptr += 5;
			if( mode == 0 ) {
				/*---------------------------------------------------------*
        * encode 20ms ACELP frame                                 *
        * acelp bits: 2+(9+6+9+6)+4+(4xICB_NBITS])+(4x7)          *
        *---------------------------------------------------------*/
				j = 7;
				/* mean energy : 2 bits */
				int2bin( prm[j], 2, ptr );
				ptr += 2;
				j++;
				for( sfr = 0; sfr < 4; sfr++ ) {
					if( ( sfr == 0 ) || ( sfr == 2 ) ) {
						n = 9;
					}
					else {
						n = 6;
					}
					/* AMR-WB closed-loop pitch lag */
					int2bin( prm[j], n, ptr );
					ptr += n;
					j++;
					int2bin( prm[j], 1, ptr );
					ptr += 1;
					j++;
					if( codec_mode == MODE_9k6 ) {
						/* 20 bits AMR-WB codebook is used */
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
					}
					else if( codec_mode == MODE_11k2 ) {
						/* 28 bits AMR-WB codebook is used */
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
						int2bin( prm[j], 5, ptr );
						ptr += 5;
						j++;
					}
					else if( codec_mode == MODE_12k8 ) {
						/* 36 bits AMR-WB codebook is used */
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
					}
					else if( codec_mode == MODE_14k4 ) {
						/* 44 bits AMR-WB codebook is used */
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
						int2bin( prm[j], 9, ptr );
						ptr += 9;
						j++;
					}
					else if( codec_mode == MODE_16k ) {
						/* 52 bits AMR-WB codebook is used */
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
						int2bin( prm[j], 13, ptr );
						ptr += 13;
						j++;
					}
					else if( codec_mode == MODE_18k4 ) {
						/* 64 bits AMR-WB codebook is used */
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
					}
					else if( codec_mode == MODE_20k ) {
						/* 72 bits AMR-WB codebook is used */
						int2bin( prm[j], 10, ptr );
						ptr += 10;
						j++;
						int2bin( prm[j], 10, ptr );
						ptr += 10;
						j++;
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 2, ptr );
						ptr += 2;
						j++;
						int2bin( prm[j], 10, ptr );
						ptr += 10;
						j++;
						int2bin( prm[j], 10, ptr );
						ptr += 10;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
						int2bin( prm[j], 14, ptr );
						ptr += 14;
						j++;
					}
					else if( codec_mode == MODE_23k2 ) {
						/* 88 bits AMR-WB codebook is used */
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
						int2bin( prm[j], 11, ptr );
						ptr += 11;
						j++;
					}
					/* AMR-WB 7 bits gains codebook */
					int2bin( prm[j], 7, ptr );
					ptr += 7;
					j++;
				}
			}    /* end of mode 0 */
			else /* mode 1 */
			{
				/* encode 20ms TCX */
				n_pack = 1;
				nbits_AVQ[0] = nbits - 56;
				AVQ_encmux( n_pack, prm + 9, prm_AVQ, nbits_AVQ, 288 / 8 );
				int2bin( prm[7], 3, ptr );
				ptr += 3;
				int2bin( prm[8], 7, ptr );
				ptr += 7;
				unpack4bits( nbits_AVQ[0], prm_AVQ, ptr );
			} /* end of mode 1 */
			k++;
		} /* end of mode 0/1 */
		else if( mode == 2 ) {
			/* encode and multiplex 40ms TCX */
			n_pack = 2;
			nbits_AVQ[0] = nbits - 26;
			nbits_AVQ[1] = nbits - 30 - 6;
			AVQ_encmux( n_pack, prm + 9, prm_AVQ, nbits_AVQ, 576 / 8 );
			/* encode first 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2; /* +2 because of mode */
			int2bin( prm[0], 8, ptr );
			ptr += 8;
			int2bin( prm[1], 8, ptr );
			ptr += 8;
			int2bin( prm[7], 3, ptr );
			ptr += 3;
			int2bin( prm[8], 7, ptr );
			ptr += 7;
			j = unpack4bits( nbits_AVQ[0], prm_AVQ, ptr );
			k++;
			/* encode second 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2;
			int2bin( prm[2], 6, ptr );
			ptr += 6;
			int2bin( prm[3], 7, ptr );
			ptr += 7;
			int2bin( prm[4], 7, ptr );
			ptr += 7;
			int2bin( prm[5], 5, ptr );
			ptr += 5;
			int2bin( prm[6], 5, ptr );
			ptr += 5;
			int2bin( prm[8] >> 1, 6, ptr );
			ptr += 6; /* 6 bits redundancy for TCX gain */
			unpack4bits( nbits_AVQ[1], prm_AVQ + j, ptr );
			k++;
		} /* end of mode 2 */
		else if( mode == 3 ) {
			/* encode and multiplex 80ms TCX */
			n_pack = 4;
			nbits_AVQ[0] = nbits - 23;
			nbits_AVQ[1] = nbits - 9 - 3;
			nbits_AVQ[2] = nbits - 12 - 3;
			nbits_AVQ[3] = nbits - 12 - 3;
			AVQ_encmux( n_pack, prm + 9, prm_AVQ, nbits_AVQ, 1152 / 8 );
			/* encode first 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2; /* +2 because of mode */
			int2bin( prm[0], 8, ptr );
			ptr += 8;
			int2bin( prm[1], 8, ptr );
			ptr += 8;
			int2bin( prm[8], 7, ptr );
			ptr += 7;
			j = unpack4bits( nbits_AVQ[0], prm_AVQ, ptr );
			k++;
			/* encode second 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2;
			int2bin( prm[2], 6, ptr );
			ptr += 6;
			int2bin( prm[7], 3, ptr );
			ptr += 3;
			/* write 3 parity check bits */
			bit = ( ( prm[8] >> 6 ) & 0x01 ) ^ ( ( prm[8] >> 3 ) & 0x01 );
			parity = bit << 2;
			bit = ( ( prm[8] >> 5 ) & 0x01 ) ^ ( ( prm[8] >> 2 ) & 0x01 );
			parity += bit << 1;
			bit = ( ( prm[8] >> 4 ) & 0x01 ) ^ ( ( prm[8] >> 1 ) & 0x01 );
			parity += bit;
			int2bin( parity, 3, ptr );
			ptr += 3;
			j += unpack4bits( nbits_AVQ[1], prm_AVQ + j, ptr );
			k++;
			/* encode third 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2;
			int2bin( prm[3], 7, ptr );
			ptr += 7;
			int2bin( prm[5], 5, ptr );
			ptr += 5;
			/* 3 bits of the TCX gain */
			int2bin( prm[8] >> 4, 3, ptr );
			ptr += 3;
			j += unpack4bits( nbits_AVQ[2], prm_AVQ + j, ptr );
			k++;
			/* encode fourth 20ms packet */
			ptr = serial + ( k * nbits_pack ) + 2;
			int2bin( prm[4], 7, ptr );
			ptr += 7;
			int2bin( prm[6], 5, ptr );
			ptr += 5;
			/* 3 bits of the TCX gain */
			int2bin( ( prm[8] >> 1 ) & 0x07, 3, ptr );
			ptr += 3;
			unpack4bits( nbits_AVQ[3], prm_AVQ + j, ptr );
			k++;
		} /* end of mode 3 */
	}     /* end of while k < NB_DIV */
	return;
}
void enc_prm_hf(
    int   mod[],     /* (i) : frame mode (mode[4], 4 division) */
    int   param[],   /* (i) : parameters                       */
    short serial[],  /* (o) : serial bits stream               */
    int   nbits_pack /* (i) : number of bits per packet of 20ms*/
)
{
	int    i, k, mode, nbits, *prm;
	short *ptr;
	/* bits per 20ms packet */
	nbits = NBITS_BWE / 4;
	k = 0;
	while( k < NB_DIV ) {
		mode = mod[k];
		/* set pointer to parameters */
		prm = param + ( k * NPRM_BWE_DIV );
		/* encode first 20ms packet */
		ptr = serial + ( ( k + 1 ) * nbits_pack ) - nbits;
		int2bin( prm[0], 2, ptr );
		ptr += 2;
		int2bin( prm[1], 7, ptr );
		ptr += 7;
		int2bin( prm[2], 7, ptr );
		ptr += 7;
		k++;
		if( mode == 2 ) {
			/* encode second 20ms packet */
			ptr = serial + ( ( k + 1 ) * nbits_pack ) - nbits;
			for( i = 3; i <= 10; i++ ) {
				int2bin( prm[i], 2, ptr );
				ptr += 2;
			}
			k++;
		}
		else if( mode == 3 ) {
			/* encode second 20ms packet */
			ptr = serial + ( ( k + 1 ) * nbits_pack ) - nbits;
			for( i = 3; i <= 10; i++ ) {
				int2bin( prm[i] >> 1, 2, ptr );
				ptr += 2;
			}
			k++;
			/* encode third 20ms packet */
			ptr = serial + ( ( k + 1 ) * nbits_pack ) - nbits;
			for( i = 11; i <= 18; i++ ) {
				int2bin( prm[i] >> 1, 2, ptr );
				ptr += 2;
			}
			k++;
			/* encode fourth 20ms packet */
			ptr = serial + ( ( k + 1 ) * nbits_pack ) - nbits;
			for( i = 3; i <= 18; i++ ) {
				int2bin( prm[i], 1, ptr );
				ptr += 1;
			}
			k++;
		}
	} /* end while (k < NB_DIV) */
	return;
}
static int unpack4bits( int nbits, int *prm, short *ptr )
{
	int i;
	i = 0;
	while( nbits > 4 ) {
		int2bin( prm[i], 4, ptr );
		ptr += 4;
		nbits -= 4;
		i++;
	}
	int2bin( prm[i], nbits, ptr );
	i++;
	return ( i );
}
