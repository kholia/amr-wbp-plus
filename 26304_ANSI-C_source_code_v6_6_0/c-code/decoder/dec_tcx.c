#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* local function */

static int find_mpitch( float xri[], int lg );
void       decoder_tcx(
          int                 prm[],       /* input:  parameters              */
          int                 nbits_AVQ[], /* input:  nbits in parameters of AVQ */
          float               A[],         /* input:  coefficients NxAz[M+1]  */
          int                 L_frame,     /* input:  frame length            */
          int                 bad_frame[], /* input:  Bad frame indicator     */
          float               exc[],       /* output: exc[-lg..lg]            */
          float               synth[],     /* in/out: synth[-M..lg]           */
          Decoder_State_Plus *st )         /* i/o : coder memory state        */
{
	int    i, k, i_subfr, index, lg, lext, bfi, n_pack;
	float  tmp, gain, fac_ns;
	float *p_A, Ap[M + 1];
	float *xri, *xnq;
	float  window[128 + 128];
	int    any_loss, num_loss;
	xri = exc;
	xnq = synth;
	if( L_frame == ( L_FRAME_PLUS / 4 ) ) {
		n_pack = 1;
	}
	else if( L_frame == ( L_FRAME_PLUS / 2 ) ) {
		n_pack = 2;
	}
	else if( L_frame == ( L_FRAME_PLUS ) ) {
		n_pack = 4;
	}
	else {
		printf( "decoder_tcx: Invalid L_frame, exiting\n" );
		exit( 1 );
	}
	any_loss = 0;
	num_loss = 0;
	for( i = 0; i < n_pack; i++ ) {
		any_loss |= bad_frame[i];
		num_loss += bad_frame[i];
	}
	/*------ set length of overlap (lext) and length of encoded frame (lg) -----*/
	lext = L_OVLP;
	if( L_frame == ( L_FRAME_PLUS / 2 ) ) {
		lext = L_OVLP / 2;
	}
	else if( L_frame == ( L_FRAME_PLUS / 4 ) ) {
		lext = L_OVLP / 4;
	}
	lg = L_frame + lext;
	/*------ initialize window of TCX target : two halves of a squared
    root of a hanning window the first part is as long as the overlap
    length of past frame (ovlp_size) the second part is as long as
    lext -----*/
	/*----------- DECODE TCX PARAMETERS AND SYNTHESIZE LOWER-BAND SIGNAL (synth[]) ----------*/
	/* if we lost all packets (i.e. 1 packet of TCX-20 ms, 2 packets of
     the TCX-40 ms or 4 packets of the TCX-80ms), we lost the whole
     coded frame extrapolation strategy: repeat lost excitation and
     use extrapolated ISFs */
	if( ( n_pack == 1 ) && ( bad_frame[0] != 0 ) ) {
		/* repeat past excitation */
		k = st->pitch_tcx;
		for( i = 0; i < L_frame; i++ ) {
			exc[i] = 0.7f * exc[i - k];
		}
		mvr2r( synth - M, synth, M );
		p_A = A;
		for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR ) {
			E_UTIL_synthesis( p_A, &exc[i_subfr], &synth[i_subfr + M], L_SUBFR, &synth[i_subfr], 0 );
			E_LPC_a_weight( p_A, Ap, GAMMA1, M );
			E_UTIL_residu( Ap, &synth[i_subfr + M], &xnq[i_subfr], L_SUBFR );
			p_A += ( M + 1 );
		}
		tmp = st->mem_wsyn;
		E_UTIL_deemph( xnq, TILT_FAC, L_frame, &tmp );
		tmp = 0.7f * st->wsyn_rms;
		st->wsyn_rms = tmp;
		for( i = 0; i < L_frame; i++ ) {
			if( xnq[i] > tmp ) {
				xnq[i] = tmp;
			}
			else {
				if( xnq[i] < -tmp ) {
					xnq[i] = -tmp;
				}
			}
		}
		E_UTIL_f_preemph( xnq, TILT_FAC, L_frame, &( st->mem_wsyn ) );
		p_A = A;
		for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR ) {
			E_LPC_a_weight( p_A, Ap, GAMMA1, M );
			E_UTIL_synthesis( Ap, &xnq[i_subfr], &synth[i_subfr], L_SUBFR, &synth[i_subfr - M], 0 );
			p_A += ( M + 1 );
		}
		/* zero windowed overlap (wovlp) for next frame */
		for( i = 0; i < lext; i++ ) {
			st->wovlp[i] = 0;
		}
	}
	else {
		/*----- synthesize windowed TCX target (inverse FFT of decoded
      spectrum + noise fill-in) using available info  ------*/
		/* decode noise level (fac_ns) (stored in 2nd packet on TCX80) */
		index = *prm++;
		if( n_pack == 4 ) {
			bfi = bad_frame[1];
		}
		else {
			bfi = bad_frame[0];
		}
		if( bfi != 0 ) {
			index = 0;
		}
		fac_ns = 0.1f * ( 8.0f - ( (float)index ) );
		/* read index of global TCX gain : erased bits are set to 1 to
       minimize global gain if bad_frame set erased bits to default
       value (0 or 1 depending on parameter) */
		index = *prm++;
		/* decode parameters of multi-rate lattice VQ */
		AVQ_demuxdec( n_pack, prm, nbits_AVQ, xri, lg / 8, bad_frame );
		if( any_loss ) {
			adapt_low_freq_deemph_ecu( xri, lg, st );
			reconst_spect( xri, st->old_xri, n_pack, bad_frame, lg, st->last_mode, xnq );
			/* generate random excitation buffer */
			set_zero( xnq, lg );
			rnd_ph16( &( st->seed_tcx ), &xnq[lg / 6], lg - ( lg / 6 ) );
			/*----------------------------------------------*
       * noise fill-in on unquantized subvector       *
       * injected only from 1066Hz to 6400Hz.         *
       *----------------------------------------------*/
			for( k = lg / 6; k < lg; k += 8 ) {
				tmp = 0.0;
				for( i = k; i < k + 8; i++ ) {
					tmp += xri[i] * xri[i];
				}
				if( tmp == 0.0 ) {
					for( i = k; i < k + 8; i++ ) {
						xri[i] = fac_ns * xnq[i];
					}
				}
			}
		}
		else {
			/* generate random excitation buffer */
			set_zero( xnq, lg );
			rnd_ph16( &( st->seed_tcx ), &xnq[lg / 6], lg - ( lg / 6 ) );

			/*----------------------------------------------*
       * noise fill-in on unquantized subvector       *
       * injected only from 1066Hz to 6400Hz.         *
       *----------------------------------------------*/
			for( k = lg / 6; k < lg; k += 8 ) {
				tmp = 0.0;
				for( i = k; i < k + 8; i++ ) {
					tmp += xri[i] * xri[i];
				}
				if( tmp == 0.0 ) {
					for( i = k; i < k + 8; i++ ) {
						xri[i] = fac_ns * xnq[i];
					}
				}
			}
			adap_low_freq_deemph( xri, lg );
		}
		/* save tcx vector for next frame */
		mvr2r( xri, st->old_xri, lg );
		/* find pitch for bfi case */
		st->pitch_tcx = find_mpitch( xri, lg );
		/*-----------------------------------------------------------*
    * Compute inverse FFT for obtaining xnq[] without noise.    *
    * Coefficients (xri[]) order are                            *
    *    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
    * Note that last FFT element (re[n/2]) is zeroed.           *
    *-----------------------------------------------------------*/
		xri[1] = 0.0; /* freq bin at 6400 Hz zeroed */
		ifft9( xri, xnq, (short)lg );
		/*-----------------------------------------------------------*
    * decode TCX global gain, multiply xnq[] by gain.           *
    * windowing of xnq[] for TCX overlap.                       *
    *-----------------------------------------------------------*/
		bfi = 0; /* We never have a missing gain!!! */
		gain = d_gain_tcx( index, xnq, lg, bfi, &( st->wsyn_rms ) );
#ifdef GAIN_CONCEAL_OLA
		gain_conceal_ola( &gain, xnq, st->wovlp, st->ovlp_size );
#endif
		for( i = 0; i < lg; i++ ) {
			xnq[i] *= gain;
		}
		/*------ initialize window of TCX target : two halves of a squared
      root of a hanning window the first part is as long as the overlap
      length of past frame (ovlp_size) the second part is as long as
      lext -----*/
		cos_window( window, st->ovlp_size, lext );
		/* adaptive windowing on overlap (beginning and end of frame) */
		for( i = 0; i < st->ovlp_size; i++ ) {
			xnq[i] *= window[i];
		}
		for( i = 0; i < lext; i++ ) {
			xnq[i + L_frame] *= window[st->ovlp_size + i];
		}
		/*-----------------------------------------------------------*
    * TCX overlap and add.  Update memory for next overlap.     *
    *-----------------------------------------------------------*/
		for( i = 0; i < L_OVLP; i++ ) {
			xnq[i] += st->wovlp[i];
		}
		/* save overlap for next frame */
		for( i = 0; i < lext; i++ ) {
			st->wovlp[i] = xnq[i + L_frame];
		}
		for( i = lext; i < L_OVLP; i++ ) {
			st->wovlp[i] = 0.0;
		}
		st->ovlp_size = lext;
		/*-----------------------------------------------------------*
    * find excitation and synthesis                             *
    *-----------------------------------------------------------*/
		E_UTIL_f_preemph( xnq, TILT_FAC, L_frame, &( st->mem_wsyn ) );
		p_A = A;
		for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR ) {
			E_LPC_a_weight( p_A, Ap, GAMMA1, M );
			E_UTIL_synthesis( Ap, &xnq[i_subfr], &synth[i_subfr], L_SUBFR, &synth[i_subfr - M], 0 );
			E_UTIL_residu( p_A, &synth[i_subfr], &exc[i_subfr], L_SUBFR );
			p_A += ( M + 1 );
		}
	}
	return;
}
static int find_mpitch( float xri[], int lg )
{
	float tmp, max, pitch, mpitch;
	int   i, n;
	max = 0.0;
	n = 2;
	/* find maximum below 400Hz */
	for( i = 2; i < lg / 16; i += 2 ) {
		tmp = xri[i] * xri[i] + xri[i + 1] * xri[i + 1];
		if( tmp > max ) {
			max = tmp;
			n = i;
		}
	}
	pitch = (float)lg * ( 2.0f / (float)n );
	/* find pitch multiple under 20ms */
	if( pitch >= 256.0f ) {
		n = 256;
	}
	else {
		mpitch = pitch;
		while( mpitch < 256.0f ) {
			mpitch += pitch;
		}
		n = (int)floor( mpitch - pitch );
	}
	return ( n );
}
