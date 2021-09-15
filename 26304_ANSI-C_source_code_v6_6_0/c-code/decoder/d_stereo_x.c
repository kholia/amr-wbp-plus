#include "../include/amr_plus.h"
#include "../include/cod_hi_stereo.h"
#include "../include/cod_tcx_stereo.h"
#include "../include/s_util.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
* Funtion  init_decoder_stereo                                    *
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                    *
*   ->Initialization of variables for the stereo decoder.         *
*-----------------------------------------------------------------*/
void init_decoder_stereo_x( Decoder_State_Plus *st )
{
	/* Allocate memory for stereo states */
	set_zero( st->my_old_synth_2k, L_FDEL_2k + D_STEREO_TCX + 2 * ( D_NC * 5 ) / 32 );
	set_zero( st->my_old_synth_hi, 2 * D_NC );
	set_zero( st->my_old_synth, 2 * L_FDEL + 20 );
	set_zero( st->mem_left_2k, 2 * L_FDEL_2k );
	set_zero( st->mem_right_2k, 2 * L_FDEL_2k );
	set_zero( st->mem_left_hi, L_FDEL );
	set_zero( st->mem_right_hi, L_FDEL );
	set_zero( st->left.mem_d_tcx, D_NC + ( D_STEREO_TCX * 32 / 5 ) );
	set_zero( st->right.mem_d_tcx, D_NC + ( D_STEREO_TCX * 32 / 5 ) );
	init_dec_hi_stereo( st );
	set_zero( st->right.mem_d_nonc, D_NC );
	set_zero( st->left.mem_d_nonc, D_NC );
	init_tcx_stereo_decoder( st );
	st->last_stereo_mode = 0;
	st->side_rms = 0.0;
	return;
}
/*-----------------------------------------------------------------*
* Funtion decoder_stereo                                          *
* ~~~~~~~~~~~~~~~~~~~~~~                                          *
*   ->Principle stereo decoder routine (working at fs=12.8kHz).   *
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
#define DELAY_MAX_2 ( D_NC + ( D_STEREO_TCX * 32 / 5 ) )
static void delay( float signal[], int lg, int delay, float mem[] )
{
	int   i;
	float buf[DELAY_MAX_2];
	for( i = 0; i < delay; i++ ) {
		buf[i] = mem[i];
	}
	for( i = 0; i < delay; i++ ) {
		mem[i] = signal[lg - delay + i];
	}
	for( i = lg - 1; i >= delay; i-- ) {
		signal[i] = signal[i - delay];
	}
	for( i = 0; i < delay; i++ ) {
		signal[i] = buf[i];
	}
	return;
}
void decoder_stereo_x( int param[],     /* (i)	:	Codebooks indices		*/
    int                    bad_frame[], /* (i)	:	Bad frame index			*/
    float                  sig_left[],  /* (o)	:	Decoded left channel	*/
    float                  synth[],     /* (o)	:	Decoded right channel	*/
    float                  AqLF[],
    int                    StbrMode,
    int                    fscale,
    Decoder_State_Plus *   st ) /* (i/o):	Decoder states			*/
{
	float *my_synth_buf = sig_left;
	float *my_new_synth = my_synth_buf + 2 * L_FDEL;
	float  my_old_synth_2k[L_FRAME_2k + D_STEREO_TCX + L_FDEL_2k + 2 * ( D_NC * 5 ) / 32];
	float *my_new_synth_2k = my_old_synth_2k + D_STEREO_TCX + L_FDEL_2k + 2 * ( D_NC * 5 ) / 32;
	float *my_old_synth_hi = synth;
	float *my_new_synth_hi = my_old_synth_hi + 2 * D_NC;
	float *my_synth_hi_t0 = my_old_synth_hi + 2 * D_NC + L_BSP + D_BPF;
	float  old_right_2k[L_FRAME_2k + 2 * L_FDEL_2k];
	float *new_right_2k = old_right_2k + 2 * L_FDEL_2k;
	float  old_left_2k[L_FRAME_2k + 2 * L_FDEL_2k];
	float *new_left_2k = old_left_2k + 2 * L_FDEL_2k;
	/* set buffers */
	mvr2r( st->my_old_synth, my_synth_buf, 2 * L_FDEL + 20 );
	mvr2r( synth, my_new_synth + 20, L_FRAME_PLUS );
	mvr2r( st->my_old_synth_2k, my_old_synth_2k, D_STEREO_TCX + L_FDEL_2k + 2 * ( D_NC * 5 ) / 32 );
	mvr2r( st->my_old_synth_hi, my_old_synth_hi, 2 * D_NC );
	mvr2r( st->mem_left_2k, old_left_2k, 2 * L_FDEL_2k );
	mvr2r( st->mem_right_2k, old_right_2k, 2 * L_FDEL_2k );
	/* mono synth band-split*/
	/* do the lo,hi band-splitting on the mono signal */
	band_split_taligned_2k( my_new_synth, my_new_synth_2k, my_new_synth_hi, L_FRAME_PLUS );

	/* Low band */
	if( StbrMode < 0 ) {
		int i;

		for( i = 0; i < L_OVLP_2k; i++ ) {
			new_left_2k[i] = my_old_synth_2k[i] + st->mem_stereo_ovlp[i];
			new_right_2k[i] = my_old_synth_2k[i] - st->mem_stereo_ovlp[i];
		}
		for( i = L_OVLP_2k; i < L_FRAME_2k; i++ ) {
			new_left_2k[i] = new_right_2k[i] = my_old_synth_2k[i];
		}

		init_tcx_stereo_decoder( st );
		st->mem_stereo_ovlp_size = L_OVLP_2k;
		st->last_stereo_mode = 0;
		st->side_rms = 0.0;
	}
	else
		dec_tcx_stereo( my_old_synth_2k, new_left_2k, new_right_2k, param + NPRM_STEREO_HI_X * NB_DIV, bad_frame, st );

	/* High band */
	{
		float  old_right_hi[L_FRAME_PLUS + L_FDEL + HI_FILT_ORDER];
		float *new_right_hi = old_right_hi + L_FDEL;
		float  old_left_hi[L_FRAME_PLUS + L_FDEL];
		float *new_left_hi = old_left_hi + L_FDEL;
		if( StereoNbits[StbrMode] - 4 > 300 ) {
			st->filt_hi_pmsvq = &filt_hi_pmsvq7;
			st->gain_hi_pmsvq = &gain_hi_pmsvq5;
		}
		else {
			st->filt_hi_pmsvq = &filt_hi_pmsvq4;
			st->gain_hi_pmsvq = &gain_hi_pmsvq2;
		}
		mvr2r( st->mem_left_hi, old_left_hi, L_FDEL );
		mvr2r( st->mem_right_hi, old_right_hi, L_FDEL );

		/* decode the high band */
		if( StbrMode < 0 ) {
			int i;

			for( i = 0; i < L_FRAME_PLUS; i++ ) {
				new_left_hi[i] = new_right_hi[i] = my_synth_hi_t0[i - L_BSP - D_BPF];
			}

			init_dec_hi_stereo( st );
		}
		else
			dec_hi_stereo( my_synth_hi_t0, new_right_hi, new_left_hi, AqLF, param, bad_frame, fscale, st );

		delay( new_right_hi, L_FRAME_PLUS, D_NC, st->right.mem_d_nonc );
		delay( new_left_hi, L_FRAME_PLUS, D_NC, st->left.mem_d_nonc );
		/*left_hi and right_hi are time aligned here */
		/* synthesis is delayed, so delay the left_hi and right_hi */
		delay( new_right_hi, L_FRAME_PLUS, D_NC + ( D_STEREO_TCX * 32 / 5 ), st->right.mem_d_tcx );
		delay( new_left_hi, L_FRAME_PLUS, D_NC + ( D_STEREO_TCX * 32 / 5 ), st->left.mem_d_tcx );
		/* the whole stereo is delayed by D_STEREO_TCX */
		/* Join low and high frequency band*/
		/* Left channel */
		mvr2r( &my_synth_buf[L_FRAME_PLUS], st->my_old_synth, 2 * L_FDEL + 20 );
		band_join_2k( sig_left, new_left_2k, new_left_hi, L_FRAME_PLUS );
		/* Right channel */
		mvr2r( &my_old_synth_hi[L_FRAME_PLUS], st->my_old_synth_hi, 2 * D_NC );
		band_join_2k( synth, new_right_2k, new_right_hi, L_FRAME_PLUS );
		mvr2r( &old_left_hi[L_FRAME_PLUS], st->mem_left_hi, L_FDEL );
		mvr2r( &old_right_hi[L_FRAME_PLUS], st->mem_right_hi, L_FDEL );
	}
	/* Update buffers*/
	mvr2r( &my_old_synth_2k[L_FRAME_2k], st->my_old_synth_2k, L_FDEL_2k + D_STEREO_TCX + 2 * ( D_NC * 5 ) / 32 );
	mvr2r( &old_left_2k[L_FRAME_2k], st->mem_left_2k, 2 * L_FDEL_2k );
	mvr2r( &old_right_2k[L_FRAME_2k], st->mem_right_2k, 2 * L_FDEL_2k );
	return;
}
static int pack4bits_d( int nbits, short *ptr, int *prm )
{
	int i;
	i = 0;
	while( nbits > 4 ) {
		prm[i] = bin2int( 4, ptr );
		ptr += 4;
		nbits -= 4;
		i++;
	}
	prm[i] = bin2int( nbits, ptr );
	i++;
	return ( i );
}
void dec_prm_stereo_x( int bad_frame[], /* (i) : bfi for 4 frames (bad_frame[4])  */
    short                  serial[],    /* (i) : serial bits stream               */
    int                    nbits_pack,  /* (i) : number of bits per packet of 20ms*/
    int                    nbits_bwe,   /* (i) : number of BWE bits per 20ms */
    int                    param[],     /* (o) : decoded parameters               */
    int                    brMode,
    Decoder_State_Plus *   st )
{
	int    nbits, *prm;
	int    k, mod_buf[1 + NB_DIV];
	int *  mod;
	int    nbits_AVQ[NB_DIV];
	int    prm_AVQ[NBITS_MAX + N_PACK_MAX];
	short *ptr;
	float  tmp_prm[NBITS_MAX + N_PACK_MAX];
	int    i, j, n_pack;
	int    hf_bits = 6;
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
	mod = mod_buf + 1;
	mod[-1] = st->last_stereo_mode; /* previous mode */
	/*----------------------------------------------------------*
	* decode the high band parameters							*
	*----------------------------------------------------------*/
	for( k = 0; k < NB_DIV; k++ ) {
		prm = param + k * NPRM_STEREO_HI_X;
		ptr = serial + ( k + 1 ) * nbits_pack - nbits;
		if( hiband_mode == 0 ) {
			prm[0] = bin2int( 4, ptr );
			ptr += 4;
			prm[1] = bin2int( 2, ptr );
			ptr += 2;
		}
		else {
			prm[0] = bin2int( 4, ptr );
			ptr += 4;
			prm[1] = bin2int( 3, ptr );
			ptr += 3;
			prm[2] = bin2int( 5, ptr );
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
		* decode the low band parameters							*
		*----------------------------------------------------------*/
	/* decode the mode */
	for( k = 0; k < NB_DIV; k++ ) {
		ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits + 1; /* +1 reserved bit*/
		if( bad_frame[k] == 0 ) {
			mod[k] = bin2int( 2, ptr );
		}
		else {
			mod[k] = -1;
		}
	}
	/* extrapolate modes in lost packets */
	if( ( mod[0] == 3 ) || ( mod[1] == 3 ) || ( mod[2] == 3 ) || ( mod[3] == 3 ) ) {
		for( k = 0; k < NB_DIV; k++ ) {
			mod[k] = 3; /* partial loss on tcx 80 ms*/
		}
	}
	else {
		if( ( mod[0] == 2 ) || ( mod[1] == 2 ) ) {
			for( k = 0; k < 2; k++ ) {
				mod[k] = 2;
			}
		}
		if( ( mod[2] == 2 ) || ( mod[3] == 2 ) ) {
			for( k = 2; k < 4; k++ ) {
				mod[k] = 2;
			}
		}
		for( k = 0; k < NB_DIV; k++ ) {
			if( mod[k] < 0 ) {
				if( mod[k - 1] == 0 ) {
					mod[k] = 0;
				}
				else {
					mod[k] = 1;
				}
			}
		}
	}
	for( k = 0; k < NB_DIV; k++ ) {
		param[NPRM_STEREO_HI_X * NB_DIV + k] = mod[k];
	}
	k = 0;
	while( k < NB_DIV ) {
		/* set pointer to parameters */
		prm = ( param + 4 + NPRM_STEREO_HI_X * NB_DIV ) + ( k * NPRM_DIV_TCX_STEREO );
		if( ( mod[k] == 1 ) || ( mod[k] == 0 ) ) {
			n_pack = 1;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 7 - 2 - hf_bits;
			/* decode 20ms TCX */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits + 2 + 1; /* +2 for the mode +1 reserved bit*/
			prm[0] = bin2int( 7, ptr );
			ptr += 7;
			prm[1] = bin2int( 7, ptr );
			ptr += 7;
			pack4bits_d( nbits_AVQ[0], ptr, prm_AVQ );
			/* demultiplex and decode */
			AVQ_demuxdec( n_pack, prm_AVQ, nbits_AVQ, tmp_prm, TOT_PRM_20 / 8, bad_frame + k );
			/* convert to integer */
			for( i = 0; i < TOT_PRM_20; i++ ) {
				prm[i + 2] = (int)tmp_prm[i];
			}
			k++;
		} /* end of mode 0/1 */
		else if( mod[k] == 2 ) {
			/* decode and demultiplex a 40 ms frame */
			n_pack = 2;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - 7 - hf_bits;
			nbits_AVQ[1] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - 7 - hf_bits;
			/* decode first 20ms packet */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits + 2 + 1; /* +2 for the mode +1 reserved bit*/
			prm[0] = bin2int( 7, ptr );
			ptr += 7;
			j = pack4bits_d( nbits_AVQ[0], ptr, prm_AVQ );
			/* decode second 20ms packet */
			ptr = serial + ( k + 2 ) * nbits_pack - nbits + hf_bits + 2 + 1; /* +2 for the mode +1 reserved bit*/
			prm[1] = bin2int( 7, ptr );
			ptr += 7;
			pack4bits_d( nbits_AVQ[1], ptr, prm_AVQ + j );
			/* demultiplex and decode tcx parameters */
			AVQ_demuxdec( n_pack, prm_AVQ, nbits_AVQ, tmp_prm, TOT_PRM_40 / 8, bad_frame + k );
			/* convert to integer */
			for( i = 0; i < TOT_PRM_40; i++ ) {
				prm[i + 2] = (int)tmp_prm[i];
			}
			k += 2;
		} /* end of mode 2 */
		else if( mod[k] == 3 ) {
			/* encode and multiplex 80ms TCX */
			n_pack = 4;
			nbits_AVQ[0] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 2 - hf_bits;
			nbits_AVQ[1] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - hf_bits;
			nbits_AVQ[2] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 7 - 2 - hf_bits;
			nbits_AVQ[3] = ( ( StereoNbits[brMode] - 4 ) / 4 ) - 2 - hf_bits;
			/* set pointer to bit stream */
			ptr = serial + ( k + 1 ) * nbits_pack - nbits + hf_bits + 2 + 1;
			/* decode first 20 ms frame */
			prm[0] = bin2int( 7, ptr );
			ptr += 7;
			j = pack4bits_d( nbits_AVQ[0], ptr, prm_AVQ );
			/* set pointer to bit stream */
			ptr = serial + ( k + 2 ) * nbits_pack - nbits + hf_bits + 2 + 1;
			/* decode second 20 ms frame */
			j += pack4bits_d( nbits_AVQ[1], ptr, prm_AVQ + j );
			/* set pointer to bit stream */
			ptr = serial + ( k + 3 ) * nbits_pack - nbits + hf_bits + 2 + 1;
			/* decode third 20 ms frame */
			prm[1] = bin2int( 7, ptr );
			ptr += 7;
			j += pack4bits_d( nbits_AVQ[2], ptr, prm_AVQ + j );
			/* set pointer to bit stream */
			ptr = serial + ( k + 4 ) * nbits_pack - nbits + hf_bits + 2 + 1;
			/* decode forth 20 ms frame */
			j += pack4bits_d( nbits_AVQ[3], ptr, prm_AVQ + j );
			/* demultiplex and decode tcx parameters */
			AVQ_demuxdec( n_pack, prm_AVQ, nbits_AVQ, tmp_prm, TOT_PRM_80 / 8, bad_frame );
			/* convert to integer */
			for( i = 0; i < TOT_PRM_80; i++ ) {
				prm[i + 2] = (int)tmp_prm[i];
			}
			k += 4;
		}
	}
	return;
}
