#include "../include/amr_plus.h"
#include "../include/cod_hi_stereo.h"
#include "../include/cod_tcx_stereo.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
* Funtion  init_coder_stereo                                      *
* ~~~~~~~~~~~~~~~~~~~~~~~~~~                                      *
*   ->Initialization of variables for the stereo coder.           *
*-----------------------------------------------------------------*/
void init_coder_stereo_x( Coder_State_Plus *st )
{
	init_cod_hi_stereo( st );
	init_tcx_stereo_encoder( st );
	return;
}
/*-----------------------------------------------------------------*
* Funtion coder_stereo                                            *
* ~~~~~~~~~~~~~~~~~~~~                                            *
*   ->Principle stereo coder routine (working at fs=12.8kHz).     *
*                                                                 *
* Note: HF band are encoded twice (2 channels) using 0.8kbps BWE. *
*       Usage of 2xBWE for stereo provide better time domain      *
*       stereo definition in HF without increasing the bit-rate.  *
*       Another advantage is that the stereo decoder is limited   *
*       to the lower band (fs=12.8kHz) and this reduce the        *
*       overall complexity of the AMR-WB+ codec.  Also, this      *
*       solution is not dependent of the AMR-WB+ mode where many  *
*       different sampling frequencies are used (16, 24, 32 kHz). *
*-----------------------------------------------------------------*/
void coder_stereo_x( float speech_hi[], /* (i)	: Mixed channel, hi	*/
    float                  chan_hi[],   /* (i)	: Right channel, hi	*/
    float                  speech_2k[], /* (i)	: Mixed channel, lo */
    float                  chan_2k[],   /* (i)	: Right channel, lo	*/
    float                  AqLF[],
    int                    brMode,
    int                    param[], /* (o)	: Encoded parameters*/
    int                    fscale,
    Coder_State_Plus *     st ) /* (i/o): Encoder states	*/
{
	int mod[4];
	int i;
	/* param */
	/*-----------------------------------------------------------------*
	* NPRM_STEREO_HI_X*NB_DIV + 4 (mode) + Rest (TCX)				   *
	*-----------------------------------------------------------------*/
	if( StereoNbits[brMode] - 4 > 300 ) {
		st->filt_hi_pmsvq = &filt_hi_pmsvq7;
		st->gain_hi_pmsvq = &gain_hi_pmsvq5;
	}
	else {
		st->filt_hi_pmsvq = &filt_hi_pmsvq4;
		st->gain_hi_pmsvq = &gain_hi_pmsvq2;
	}
	/* encode the high band */
	cod_hi_stereo( speech_hi, chan_hi, AqLF, param, st ); /* 64 echantillon look-ahead */
	/* encode the low band */
	cod_tcx_stereo( speech_2k, chan_2k, param + 4 + NPRM_STEREO_HI_X * NB_DIV, brMode, mod, fscale, st );
	/* transmitt the mode in the parameters buffer*/
	for( i = 0; i < 4; i++ ) {
		param[i + NPRM_STEREO_HI_X * NB_DIV] = mod[i];
	}
	return;
}
static int unpack4bits_d( int nbits, int *prm, short *ptr )
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
void enc_prm_stereo_x( int param[],    /* (i) : parameters                       */
    short                  serial[],   /* (o) : serial bits stream               */
    int                    nbits_pack, /* (i) : number of bits per packet of 20ms*/
    int                    nbits_bwe,  /* (i) : number of BWE bits per 20ms  */
    int                    brMode )
{
	int    nbits, *prm;
	int    k, j;
	short *ptr;
	int    i;
	int    mod[NB_DIV], nbits_AVQ[NB_DIV];
	int    mode, n_pack;
	int    prm_AVQ[( NBITS_MAX / 4 ) + N_PACK_MAX];
	int    hf_bits = 0;
	int    hiband_mode;
	/*----------------------------------------------------------*
	* Set number of bits used for stereo (per packet of 20 ms) *
	* When stereo is transmitted, the bit ordering is:         *
	* serial: mode (2bits), core, stereo, 2xBWE(2x16bits)      *
	*----------------------------------------------------------*/
	nbits = ( StereoNbits[brMode] + ( 2 * nbits_bwe ) ) / 4;
	hiband_mode = 0;
	if( StereoNbits[brMode] - 4 > 300 ) {
		hiband_mode = 1;
	}
	/*----------------------------------------------------------*
	* Encode the high band parameters	 		   *
	*----------------------------------------------------------*/
	for( k = 0; k < NB_DIV; k++ ) {
		prm = param + k * NPRM_STEREO_HI_X;
		ptr = serial + ( k + 1 ) * nbits_pack - nbits;
		if( hiband_mode == 0 ) {
			int2bin( prm[0], 4, ptr );
			ptr += 4;
			int2bin( prm[1], 2, ptr );
			ptr += 2;
		}
		else {
			int2bin( prm[0], 4, ptr );
			ptr += 4;
			int2bin( prm[1], 3, ptr );
			ptr += 3;
			int2bin( prm[2], 5, ptr );
			ptr += 5;
		}
	}
	if( hiband_mode == 0 ) {
		hf_bits = 4 + 2;
	}
	else {
		hf_bits = 7 + 5;
	}
	/*----------------------------------------------------------*
		* Encode the low band parameters							*
		*----------------------------------------------------------*/
	/* fill up the mode */
	for( i = 0; i < NB_DIV; i++ ) {
		mod[i] = param[i + NPRM_STEREO_HI_X * NB_DIV];
	}
	k = 0;
	while( k < NB_DIV ) {
		mode = mod[k];
		/* set pointer to parameters */
		prm = ( param + 4 + NPRM_STEREO_HI_X * NB_DIV ) + ( k * NPRM_DIV_TCX_STEREO );
		if( ( mode == 1 ) || ( mode == 0 ) ) {
			/* encode 20ms TCX */
			n_pack = 1;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 2 - 7 - hf_bits;
			AVQ_encmux( n_pack, prm + 2, prm_AVQ, nbits_AVQ, TOT_PRM_20 / 8 );
			/* set pointer to bit stream */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits;
			/* encode the mode */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			int2bin( prm[0], 7, ptr );
			ptr += 7;
			int2bin( prm[1], 7, ptr );
			ptr += 7;
			unpack4bits_d( nbits_AVQ[0], prm_AVQ, ptr );
			ptr += nbits_AVQ[0];
			k++;
		} /* end of mode 0/1 */
		else if( mode == 2 ) {
			/* encode and multiplex 40ms TCX */
			n_pack = 2;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - 7 - hf_bits;
			nbits_AVQ[1] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - 7 - hf_bits;
			AVQ_encmux( n_pack, prm + 2, prm_AVQ, nbits_AVQ, TOT_PRM_40 / 8 );
			/* set pointer to bit stream */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits;
			/* encode first 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			int2bin( prm[0], 7, ptr );
			ptr += 7;
			j = unpack4bits_d( nbits_AVQ[0], prm_AVQ, ptr );
			/* set pointer to bit stream */
			ptr = serial + ( k + 2 ) * nbits_pack - nbits + hf_bits;
			/* encode second 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			int2bin( prm[1], 7, ptr );
			ptr += 7;
			unpack4bits_d( nbits_AVQ[1], prm_AVQ + j, ptr );
			k += 2;
		} /* end of mode 2 */
		else if( mode == 3 ) {
			/* encode and multiplex 80ms TCX */
			n_pack = 4;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 2 - hf_bits;
			nbits_AVQ[1] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - hf_bits;
			nbits_AVQ[2] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 2 - hf_bits;
			nbits_AVQ[3] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - hf_bits;
			AVQ_encmux( n_pack, prm + 2, prm_AVQ, nbits_AVQ, TOT_PRM_80 / 8 );
			/* set pointer to bit stream */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits;
			/* encode first 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			int2bin( prm[0], 7, ptr );
			ptr += 7;
			j = unpack4bits_d( nbits_AVQ[0], prm_AVQ, ptr );
			/* set pointer to bit stream */
			ptr = serial + ( k + 2 ) * nbits_pack - nbits + hf_bits;
			/* encode second 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			j += unpack4bits_d( nbits_AVQ[1], prm_AVQ + j, ptr );
			/* set pointer to bit stream */
			ptr = serial + ( k + 3 ) * nbits_pack - nbits + hf_bits;
			/* encode third 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			int2bin( prm[1], 7, ptr );
			ptr += 7;
			j += unpack4bits_d( nbits_AVQ[2], prm_AVQ + j, ptr );
			/* set pointer to bit stream */
			ptr = serial + ( k + 4 ) * nbits_pack - nbits + hf_bits;
			/* encode forth 20 ms frame */
			*ptr = 0;
			ptr += 1;
			int2bin( mode, 2, ptr );
			ptr += 2;
			unpack4bits_d( nbits_AVQ[3], prm_AVQ + j, ptr );
			k += 4;
		} /* end of mode 3 */
	}     /* end of while k < NB_DIV */
	return;
}
