#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
 * Funtion  init_decoder_lf                                        *
 * ~~~~~~~~~~~~~~~~~~~~~~~~                                        *
 *                                                                 *
 *   ->Initialization of variables for the decoder section.        *
 *                                                                 *
 *-----------------------------------------------------------------*/
void init_decoder_lf( Decoder_State_Plus *st )
{
	int i;

	/* Static vectors to zero */
	set_zero( st->old_exc, PIT_MAX_MAX + L_INTERPOL ); /* old excitation vector (>20ms) */
	set_zero( st->old_synth, M );
	set_zero( st->past_isfq, M );
	set_zero( st->wovlp, L_OVLP );
	/* Initialize the ISFs */
	mvr2r( (float *)isf_init, st->isfold, M );
	for( i = 0; i < L_MEANBUF; i++ ) {
		mvr2r( (float *)isf_init, &( st->isf_buf[i * M] ), M );
	}
	/* Initialize the ISPs */
	for( i = 0; i < M - 1; i++ ) {
		st->ispold[i] = (float)cos( 3.141592654 * (float)( i + 1 ) / (float)M );
	}
	st->ispold[M - 1] = 0.045f;
	st->prev_lpc_lost = 0;
	st->ovlp_size = 0;
	st->old_T0 = 64;
	st->old_T0_frac = 0;
	st->seed_ace = 0;
	st->mem_wsyn = 0.0;
	st->seed_tcx = 21845;
	st->wsyn_rms = 0.0;
	st->past_gpit = 0.0;
	st->past_gcode = 0.0;
	st->pitch_tcx = L_DIV;
	st->gc_threshold = 0.0f;
	return;
}
extern float isfnew_uq[NB_DIV * M];
/*-----------------------------------------------------------------*
 * Funtion decoder_lf                                              *
 * ~~~~~~~~~~~~~~~~~~                                              *
 *   ->Principle decoder routine (band 0 to 6400Hz).               *
 *                                                                 *
 *-----------------------------------------------------------------*/
void decoder_lf(
    int                 mod[],       /* (i)  : mode for each 20ms frame (mode[4]     */
    int                 param[],     /* (i)  : parameters                            */
    int                 nbits_AVQ[], /* (i)  : for each frame (nbits_AVQ[4])         */
    int                 codec_mode,  /* (i)  : AMR-WB+ mode (see cnst.h)             */
    int                 bad_frame[], /* (i)  : for each frame (bad_frame[4])         */
    float               AqLF[],      /* (o)  : decoded coefficients (AdLF[16])       */
    float               exc[],       /* (i/o): excitation[-(PIT_MAX_MAX+L_INTERPOL)..L_TCX-1]   */
    float               fsynth[],    /* (o)  : decoded synthesis                     */
    int                 pitch[],     /* (o)  : decoded pitch (pitch[16])             */
    float               pit_gain[],  /* (o)  : decoded pitch gain (pit_gain[16])     */
    int                 pit_adj,
    Decoder_State_Plus *st ) /* i/o  : decoder memory state             */
{
	/* LPC coefficients */
	float  ispnew[M];
	float  isfnew[M];
	float  Aq[( NB_SUBFR + 1 ) * ( M + 1 )]; /* A(z) for all subframes     */
	float *synth;
	/* Scalars */
	int i, k, mode, bfi, bfi_2nd_st, *prm;
	/* LTP parameters for high band */
	int   T_out[4];
	float p_out[4];
	int   len;
	float tmp, stab_fac;
	if( pit_adj == 0 ) {
		len = L_DIV;
	}
	else {
		i = ( ( ( pit_adj * PIT_MIN_12k8 ) + ( FSCALE_DENOM / 2 ) ) / FSCALE_DENOM ) - PIT_MIN_12k8;
		len = ( PIT_MAX_12k8 + ( 6 * i ) ) + L_INTERPOL;
		if( len < L_DIV ) {
			len = L_DIV;
		}
	}
	/* initialize bass postfilter */
	for( i = 0; i < NB_SUBFR; i++ ) {
		pitch[i] = L_SUBFR;
	}
	for( i = 0; i < NB_SUBFR; i++ ) {
		pit_gain[i] = 0.0;
	}
	/* Initialize pointers */
	synth = fsynth;
	mvr2r( st->old_synth, synth - M, M );
	k = 0;
	while( k < NB_DIV ) {
		mode = mod[k];
		bfi = bad_frame[k];
		bfi_2nd_st = 0;
		if( ( mode == 2 ) && ( bad_frame[k + 1] != 0 ) ) {
			bfi_2nd_st = ( 1 + 2 + 4 + 8 + 16 );
		}
		if( mode == 3 ) {
			if( bad_frame[k + 1] != 0 ) {
				bfi_2nd_st = 1;
			}
			if( bad_frame[k + 2] != 0 ) {
				bfi_2nd_st += ( 2 + 8 );
			}
			if( bad_frame[k + 3] != 0 ) {
				bfi_2nd_st += ( 4 + 16 );
			}
		}
		/* set pointer to parameters */
		prm = param + ( k * DEC_NPRM_DIV );
		/* decode ISFs and convert ISFs to cosine domain */
		dpisf_2s_46b( prm, isfnew, st->past_isfq, st->isfold, st->isf_buf, bfi, bfi_2nd_st, 1 );

		if( !bfi && bfi_2nd_st == 0 )
			st->seed_tcx = (short)( ( prm[0] * prm[1] + prm[4] * prm[5] ) * 256 );

		prm += 7;
		isf2isp( isfnew, ispnew, M );

		/* Check stability on isf : distance between old isf and current isf */
		tmp = 0.0f;
		for( i = 0; i < M - 1; i++ ) {
			tmp += ( isfnew[i] - st->isfold[i] ) * ( isfnew[i] - st->isfold[i] );
		}
		stab_fac = (float)( 1.25f - ( tmp / 400000.0f ) );
		if( stab_fac > 1.0f ) {
			stab_fac = 1.0f;
		}
		if( stab_fac < 0.0f ) {
			stab_fac = 0.0f;
		}
		/* ajust old isp[] following to a bad frame (to avoid overshoot) */
		if( ( st->prev_lpc_lost != 0 ) && ( bfi == 0 ) ) {
			for( i = 0; i < M - 1; i++ ) {
				if( ispnew[i] < st->ispold[i] ) {
					st->ispold[i] = ispnew[i];
				}
			}
		}
		/* - interpolate Ai in ISP domain (Aq) and save values for upper-band (Aq_lpc)
       - decode other parameters according to mode
       - set overlap size for next decoded frame (ovlp_size)
       - set mode for upper-band decoder (mod[])*/
		switch( mode ) {
		case 0:
		case 1:
			int_lpc_np1( st->ispold, ispnew, Aq, 4, M );
			mvr2r( Aq, &AqLF[( k * 4 ) * ( M + 1 )], 5 * ( M + 1 ) );

			if( mode == 0 ) /* 20-ms ACELP */
			{
				/* set global variables used in dec_ace() even for bfi=1*/
				decoder_acelp( prm, Aq, L_DIV, codec_mode, bad_frame[k], &exc[k * L_DIV], &synth[k * L_DIV], T_out, p_out, &pitch[k * 4], &pit_gain[k * 4], pit_adj, stab_fac, st );
				/* average integer pitch-lag for high band coding */
				st->last_mode = 0;
			}
			else { /* 20+2.5-ms TCX */
				decoder_tcx( prm, &nbits_AVQ[k], Aq, L_DIV, bad_frame + k, &exc[k * L_DIV], &synth[k * L_DIV], st );
				st->last_mode = 1;
			}
			k++;
			break;
		case 2: /* 40+5-ms TCX */
			int_lpc_np1( st->ispold, ispnew, Aq, 8, M );
			mvr2r( Aq, &AqLF[( k * 4 ) * ( M + 1 )], 9 * ( M + 1 ) );
			decoder_tcx( prm, &nbits_AVQ[k], Aq, 2 * L_DIV, bad_frame + k, &exc[k * L_DIV], &synth[k * L_DIV], st );
			st->last_mode = 2;
			k += 2;
			break;
		case 3: /* 80+10-ms TCX */
			int_lpc_np1( st->ispold, ispnew, Aq, 16, M );
			mvr2r( Aq, &AqLF[( k * 4 ) * ( M + 1 )], 17 * ( M + 1 ) );
			decoder_tcx( prm, &nbits_AVQ[k], Aq, 4 * L_DIV, bad_frame, &exc[k * L_DIV], &synth[k * L_DIV], st );
			st->last_mode = 3;
			k += 4;
			break;
		default:
			printf( "decoder error: mode > 3!\n" );
			exit( 0 );
		}
		st->prev_lpc_lost = bfi;
		/* update ispold[] and isfold[] for the next frame */
		mvr2r( ispnew, st->ispold, M );
		mvr2r( isfnew, st->isfold, M );
	}
	/*------ update signal for next superframe :
   shift exc[], synth[] and synth_hf[] to the left by L_FRAME_PLUS -----*/
	mvr2r( &synth[L_FRAME_PLUS - M], st->old_synth, M );
	/* output exc and synth */
	mvr2r( synth, fsynth, L_FRAME_PLUS );
	return;
}
