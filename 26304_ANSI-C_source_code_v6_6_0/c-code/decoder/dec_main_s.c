#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define DELAY_MAX ( D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 )
static void delay( float signal[], int lg, int delay, float mem[] )
{
	int   i;
	float buf[DELAY_MAX];
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
/*-----------------------------------------------------------------*
*      Initialization of variables for the decoder section        *
*-----------------------------------------------------------------*/
void init_decoder_amrwb_plus(
    Decoder_State_Plus *st,
    int                 num_chan,
    int                 fscale,
    short               full_reset )
{

	/* initialize memories (stereo part) */
	st->left.over_frac = 0;
	st->right.over_frac = 0;
	if( fscale == 0 ) {
		st->left.mean_isf_hf = (float *)mean_isf_hf_low_rate;
		st->left.dico1_isf_hf = (float *)dico1_isf_hf_low_rate;
		st->right.mean_isf_hf = (float *)mean_isf_hf_low_rate;
		st->right.dico1_isf_hf = (float *)dico1_isf_hf_low_rate;
	}
	else {
		st->left.mean_isf_hf = (float *)mean_isf_hf_12k8;
		st->left.dico1_isf_hf = (float *)dico1_isf_hf_12k8;
		st->right.mean_isf_hf = (float *)mean_isf_hf_12k8;
		st->right.dico1_isf_hf = (float *)dico1_isf_hf_12k8;
	}
	if( full_reset ) {
		init_bass_postfilter( st );
		/* init filters memories */
		set_zero( st->left.mem_oversamp, L_MEM_JOIN_OVER );
		set_zero( st->right.mem_oversamp, L_MEM_JOIN_OVER );
		set_zero( st->left.mem_oversamp_hf, 2 * L_FILT );
		set_zero( st->right.mem_oversamp_hf, 2 * L_FILT );
		set_zero( st->left.old_synth_hf, ( D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ) );
		set_zero( st->right.old_synth_hf, ( D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ) );
		init_decoder_hf( &( st->left ) );
		init_decoder_hf( &( st->right ) );
		set_zero( st->mem_sig_out, 4 );
		set_zero( st->left.mem_sig_out, 4 );
		set_zero( st->right.mem_sig_out, 4 );
		/* initialize memories (mono part) */
		st->last_mode = 0; /* set last mode to acelp */
		st->mem_deemph = 0.0;
		set_zero( st->h, ECU_WIEN_ORD + 1 );
		init_decoder_lf( st );
		if( num_chan == 2 ) {
			init_decoder_stereo_x( st );
		}
		set_zero( st->old_xri, L_TCX );
		set_zero( st->mem_gain_code, 4 );
		set_zero( st->mem_lpc_hf + 1, MHF );
		st->mem_lpc_hf[0] = 1.0f;
		st->mem_gain_hf = 0.0f;
		st->ramp_state = 64;
	}
	return;
}
/*-----------------------------------------------------------------*
 *                  Principal decoder routine                      *
 *-----------------------------------------------------------------*/
int decoder_amrwb_plus(                  /* output: number of sample processed */
    int                 codec_mode,      /* input: AMR-WB+ mode (see cnst.h)         */
    short               serial[],        /* input: serial parameters (4x20ms)        */
    int                 bad_frame[],     /* input: bfi (bad_frame[4])                */
    int                 L_frame,         /* input: frame size of synthesis           */
    int                 n_channel,       /* input: 1 or 2 (mono/stereo)              */
    float               channel_right[], /* (o): used on mono and stereo         */
    float               channel_left[],  /* (o): used on stereo only             */
    Decoder_State_Plus *st,              /* i/o:  decoder memory state              */
    int                 fscale,
    int                 StbrMode,
    int                 mono_dec_stereo,
    short               upsamp_fscale )
{
	float AqLF[( NB_SUBFR + 1 ) * ( M + 1 )]; /* LPC coefficients for band 0..6400Hz */
	float sig_left[L_FRAME_PLUS + 2 * L_FDEL + 20];
	float synth_buf[L_BSP + L_TCX];
	int   pitch[NB_SUBFR];
	float pit_gain[NB_SUBFR];
	int   nbits_AVQ[4];
	int   param[DEC_NPRM_DIV * NB_DIV];
	int   prm_stereo[MAX_NPRM_STEREO_DIV * NB_DIV]; /* see cnst.h */
	int   prm_hf_left[NPRM_BWE_DIV * NB_DIV];
	int   prm_hf_right[NPRM_BWE_DIV * NB_DIV];
	/* Scalars */
	int    i, k;
	int    mod_buf[1 + NB_DIV], *mod;
	int    nbits_pack;
	float *synth;
	int    bad_frame_hf[4];
	int    n_loss, any_tcx80;
	int    nb_samp, fac_fs;
	synth = synth_buf + L_FILT_2k;

	/* for 20-ms packetization, divide by 4 the 80-ms bitstream */
	nbits_pack = ( NBITS_CORE[codec_mode] + NBITS_BWE ) / 4;

	mono_dec_stereo = 0;
	if( StbrMode >= 0 ) {
		if( n_channel == 1 )
			mono_dec_stereo = 1;

		nbits_pack += ( StereoNbits[StbrMode] + NBITS_BWE ) / 4;
	}
	/*------------------------------------------------------------------------*
   * read modes (somes modes may be unreadable if the related frame is lost)*
   *------------------------------------------------------------------------*
   * mode 0 = ACELP 20ms   --> mode[] = 0,x,x,x  ..  x,x,x,0                *
   * mode 1 = TCX 20ms     --> mode[] = 1,x,x,x  ..  x,x,x,1                *
   * mode 2 = TCX 40ms     --> mode[] = 2,2,x,x  or  x,x,2,2                *
   * mode 3 = TCX 80ms     --> mode[] = 3,3,3,3                             *
   *------------------------------------------------------------------------*/
	mod = mod_buf + 1;
	mod[-1] = st->last_mode; /* previous mode */
	n_loss = 0;
	for( k = 0; k < NB_DIV; k++ ) {
		if( bad_frame[k] == 0 ) {
			mod[k] = bin2int( 2, &serial[k * nbits_pack] );
		}
		else {
			mod[k] = -1;
		}
		n_loss += bad_frame[k];
	}
	for( i = 0; i < 4; i++ ) {
		bad_frame_hf[i] = bad_frame[i];
	}
	any_tcx80 = ( mod[0] == 3 );
	for( i = 1; i < 4; i++ ) {
		any_tcx80 |= ( mod[i] == 3 );
	}
	if( ( n_loss > 2 ) && ( any_tcx80 ) ) {
		for( i = 0; i < 4; i++ ) {
			mod[i] = -1;
			bad_frame[i] = 1;
		}
	}
	/* extrapolate mode in case of packet erasures and fix bit errors */
	if( ( mod[0] == 3 ) || ( mod[1] == 3 ) || ( mod[2] == 3 ) || ( mod[3] == 3 ) ) {
		/* handle loss of one or several TCX-80 packets */
		for( k = 0; k < NB_DIV; k++ ) {
			mod[k] = 3;
		}
	}
	else {
		/* handle loss of a TCX-40 packet */
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
		/* handle loss of an ACELP or TCX-20 packet */
		for( k = 0; k < NB_DIV; k++ ) {
			if( mod[k] < 0 ) {
				/* repeat the previous mode:
	   if previous mode = ACELP        -> ACELP
	   if previous mode = TCX-20/40/80 -> TCX-20
	   notes:
	   - ACELP is not allowed after TCX (no pitch information to reuse)
	   - TCX-40 is not allowed in the mode repetition to keep the logic simple */
				if( mod[k - 1] == 0 ) {
					mod[k] = 0;
				}
				else {
					mod[k] = 1;
				}
			}
		}
	}
	/*----------- DECODE BIT-STREAM ----------*/
	dec_prm( mod, bad_frame, serial, nbits_pack, codec_mode, param, nbits_AVQ );
	if( ( n_channel == 2 ) && ( StbrMode >= 0 ) ) {
		dec_prm_stereo_x( bad_frame_hf, serial, nbits_pack, NBITS_BWE, prm_stereo, StbrMode, st );
	}
	dec_prm_hf( mod, bad_frame, serial, nbits_pack, prm_hf_right );
	if( n_channel == 2 || mono_dec_stereo == 1 ) {
		if( StbrMode < 0 )
			mvi2i( prm_hf_right, prm_hf_left, NPRM_BWE_DIV * NB_DIV );
		else
			dec_prm_hf( mod, bad_frame, serial - NBITS_BWE / 4, nbits_pack, prm_hf_left );
	}
	decoder_amrwb_plus_1( channel_right, channel_left, mod, param, prm_hf_right, prm_hf_left, nbits_AVQ, codec_mode, bad_frame, bad_frame_hf, AqLF, synth, pitch, pit_gain, st, n_channel, L_frame, fscale, mono_dec_stereo );

	if( n_channel == 2 ) {
		decoder_stereo_x( prm_stereo, bad_frame_hf, sig_left, synth, AqLF, StbrMode, fscale, st );
		if( fscale == 0 ) {
			hp50_12k8( sig_left, L_FRAME_PLUS, st->right.mem_sig_out, fscale );
			hp50_12k8( synth, L_FRAME_PLUS, st->left.mem_sig_out, fscale );
		}
		if( fscale == 0 ) {
			oversamp_12k8( sig_left, channel_left, L_frame, st->left.mem_oversamp, 0, 1 );
			oversamp_12k8( synth, channel_right, L_frame, st->right.mem_oversamp, 0, 1 );
			nb_samp = L_frame;
		}
		else {
			/* band join and oversampling (HF into channel_left/right) */
			if( L_frame == L_FRAME44k )
				fac_fs = upsamp_fscale;
			else
				fac_fs = ( L_FRAME48k * upsamp_fscale ) / L_frame;

			join_over_12k8( sig_left, channel_left, L_FRAME_PLUS, channel_left, L_frame, fac_fs, st->left.mem_oversamp, &( st->left.over_frac ) );

			nb_samp = join_over_12k8( synth, channel_right, L_FRAME_PLUS, channel_right, L_frame, fac_fs, st->right.mem_oversamp, &( st->right.over_frac ) );
		}
	}
	else {
		if( fscale == 0 ) {
			oversamp_12k8( synth, channel_right, L_frame, st->right.mem_oversamp, 0, 1 );
			nb_samp = L_frame;
		}
		else {
			/* band join and oversampling (HF into channel_left/right) */
			if( L_frame == L_FRAME44k )
				fac_fs = upsamp_fscale;
			else
				fac_fs = ( L_FRAME48k * upsamp_fscale ) / L_frame;

			nb_samp = join_over_12k8( synth, channel_right, L_FRAME_PLUS, channel_right, L_frame, fac_fs, st->right.mem_oversamp, &( st->right.over_frac ) );
		}
	}
	/* update for next superframe */
	st->last_mode = mod[NB_DIV - 1];
	return ( nb_samp );
}
void decoder_amrwb_plus_1( float *chan_right,
    float *                       chan_left,
    int *                         mod,
    int *                         param,
    int *                         prm_hf_right,
    int *                         prm_hf_left,
    int *                         nbits_AVQ,
    int                           codec_mode,
    int *                         bad_frame,
    int *                         bad_frame_hf,
    float *                       AqLF,
    float *                       synth,
    int *                         pitch,
    float *                       pit_gain,
    Decoder_State_Plus *          st,
    int                           n_channel,
    int                           L_frame,
    int                           fscale,
    int                           mono_dec_stereo )
{
	float  exc_buf[PIT_MAX_MAX + L_INTERPOL + L_TCX];
	float *exc;

	exc = exc_buf + PIT_MAX_MAX + L_INTERPOL;
	mvr2r( st->old_exc, exc_buf, PIT_MAX_MAX + L_INTERPOL );
	decoder_lf( mod, param, nbits_AVQ, codec_mode, bad_frame, AqLF, exc, synth, pitch, pit_gain, fscale, st );
	mvr2r( exc_buf + L_FRAME_PLUS, st->old_exc, PIT_MAX_MAX + L_INTERPOL );
	E_UTIL_deemph( synth, PREEMPH_FAC, L_FRAME_PLUS, &( st->mem_deemph ) );
	bass_postfilter( synth, pitch, pit_gain, synth, fscale, st );
	if( fscale == 0 ) {
		if( n_channel == 1 ) {
			/* high pass filter only in the mono case to avoid phase mismatch */
			hp50_12k8( synth, L_FRAME_PLUS, st->mem_sig_out, fscale );
		}
	}
	if( L_frame > L_FRAME8k ) {
		decoder_hf( mod, prm_hf_right, prm_hf_left, mono_dec_stereo, bad_frame_hf, AqLF, exc, chan_right, st->mem_lpc_hf, &( st->mem_gain_hf ), &( st->ramp_state ), &( st->right ) );

		if( n_channel == 1 ) {
			delay( chan_right, L_FRAME_PLUS, ( fscale == 0 ) ? ( DELAY_PF ) : ( DELAY_PF + L_SUBFR ), st->right.old_synth_hf );
		}
		else {
			delay( chan_right, L_FRAME_PLUS, ( fscale == 0 ) ? ( D_BPF + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ) : ( D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ), st->right.old_synth_hf );
		}
		if( n_channel == 2 ) {
			decoder_hf( mod, prm_hf_left, prm_hf_right, mono_dec_stereo, bad_frame_hf, AqLF, exc, chan_left, st->mem_lpc_hf, &( st->mem_gain_hf ), &( st->ramp_state ), &( st->left ) );
			delay( chan_left, L_FRAME_PLUS, ( fscale == 0 ) ? ( D_BPF + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ) : ( D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5 ), st->left.old_synth_hf );
		}
		if( fscale == 0 ) {
			oversamp_12k8( chan_right, chan_right, L_frame, st->right.mem_oversamp_hf, 1, 0 );
			if( n_channel == 2 ) {
				oversamp_12k8( chan_left, chan_left, L_frame, st->left.mem_oversamp_hf, 1, 0 );
			}
		}
	}
	else {
		set_zero( chan_right, L_frame );
		set_zero( chan_left, L_frame );
	}
	return;
}
