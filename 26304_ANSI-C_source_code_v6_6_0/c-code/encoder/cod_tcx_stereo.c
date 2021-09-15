#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* prototypes*/
void init_tcx_stereo_encoder( Coder_State_Plus *st )
{
	st->mem_stereo_ovlp_size = 0;
	set_zero( st->mem_stereo_ovlp, L_OVLP_2k );
}
int q_gain_pan( /* output: return quantization index */
    float *gain /* in/out: quantized gain            */
)
{
	int index;
	index = (int)floor( ( *gain + 2.0f ) * ( 32.0f ) + 0.5f );
	if( index < 0 )
		index = 0;
	if( index > 127 )
		index = 127;
	*gain = (float)index / 32.0f - 2.0f;
	return ( index ); /* 0...127*/
}
void ctcx_stereo( float side[], /* input:  speech[-M..lg]          */
    float               mono[],
    float               synth[],   /* in/out: synth[-M..lg]           */
    float               wovlp[],   /* i/o:    wovlp[0..127]           */
    int                 ovlp_size, /* input:  0, 64 or 128 (0=acelp)  */
    int                 L_frame,   /* input:  frame length            */
    int                 nb_bits,   /* input:  number of bits allowed  */
    int                 prm[],     /* output: tcx parameters          */
    int                 pre_echo )
{
	int   i, k, i_subfr, lg;
	int   lext = 32;
	float tmp, gain, fac_ns, gain_pan;
	float xri[L_TCX_LB];
	float xn[L_TCX_LB];
	float wm[L_TCX_LB];
	float xnq[L_TCX_LB];
#ifndef COS_FAC
	float window[L_TCX_LB];
#else
	float tmpfloat;
#endif
	float gain_shap[8];
	/*------ set length of overlap (lext) and length of encoded frame (lg) -----*/
	switch( L_frame ) {
	case 40:
		lext = 8;
		break;
	case 80:
		lext = 16;
		break;
	case 160:
		lext = 32;
		break;
	};
	lg = L_frame + lext;
	/* built window for overlaps section */
#ifndef COS_FAC
	cos_window( window, ovlp_size, lext );
#endif
	/* reduce by the correlation with the mono  */
	for( i = 0; i < lg; i++ ) {
		xn[i] = side[i];
		wm[i] = mono[i];
	}
	/* xn[] windowing for TCX overlap and correlation */
	for( i = 0; i < ovlp_size; i++ ) {
#ifndef COS_FAC
		xn[i] *= window[i];
		wm[i] *= window[i];
#else
		tmpfloat = cos_fac( i, ovlp_size, lext );
		xn[i] *= tmpfloat;
		wm[i] *= tmpfloat;
#endif
	}

	for( i = 0; i < lext; i++ ) {
#ifndef COS_FAC
		xn[L_frame + i] *= window[ovlp_size + i];
		wm[L_frame + i] *= window[ovlp_size + i];
#else
		tmpfloat = cos_fac( ovlp_size + i, ovlp_size, lext );
		xn[L_frame + i] *= tmpfloat;
		wm[L_frame + i] *= tmpfloat;
#endif
	}

	if( pre_echo ) {
		/* compensate for the gain */
		tmp = 0.0f;
		for( i = 0; i < lg; i += 16 ) {
			gain_shap[i / 16] = 0.001f;
			for( k = 0; k < 16; k++ ) {
				gain_shap[i / 16] += wm[i + k] * wm[i + k];
			}
			/* average log gain in frame */
			tmp += (float)log10( gain_shap[i / 16] );
		}
		tmp /= (float)( lg / 16 );
		for( i = 0; i < lg; i += 16 ) {
			gain_shap[i / 16] = (float)sqrt( gain_shap[i / 16] * pow( 10.0f, -tmp ) );

			gain_shap[i / 16] = my_min( 2.0, my_max( 0.5, gain_shap[i / 16] ) );

			for( k = 0; k < 16; k++ ) {
				xn[i + k] /= gain_shap[i / 16];

				wm[i + k] /= gain_shap[i / 16];
			}
		}
	}
	/* compute the optimal panning gain */
	gain_pan = get_gain( xn, wm, lg );
	/* do not amplify the mono */
	/* quantize the panning gain (would be better if predictive scalar quantizer is used)*/
	nb_bits -= 7;
	prm[0] = q_gain_pan( &gain_pan );
	/* compute the target */
	for( i = 0; i < lg; i++ ) {
		xn[i] -= gain_pan * wm[i];
	}
	/*-----------------------------------------------------------*
	* Compute the FFT of xn[].                                  *
	* Coefficients (xri[]) order are                            *
	*    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
	* Note that last FFT element (re[n/2]) is zeroed.           *
	*-----------------------------------------------------------*/
	fft3( xn, xri, (short)lg );
	xri[1] = 0.0; /* freq bin at 1000 Hz zeroed */
	/*-----------------------------------------------------------*
	* Spectral algebraic quantization                           *
	* with adaptive low frequency emphasis/deemphasis.          *
	* Noise factor is the average level of unquantized freq.    *
	*-----------------------------------------------------------*/
	adap_low_freq_emph( xri, 4 * lg );
	nb_bits -= ( 7 ); /* gain = 7 bits */
	fac_ns = AVQ_cod( xri, prm + 2, nb_bits, lg / 8 );
	for( i = 0; i < lg; i++ ) {
		xri[i] = (float)prm[i + 2];
	}
	/*-----------------------------------------------------------*
	* Compute inverse FFT for obtaining xnq[] without noise.    *
	* Coefficients (xri[]) order are                            *
	*    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
	* Note that last FFT element (re[n/2]) is zeroed.           *
	*-----------------------------------------------------------*/
	adap_low_freq_deemph( xri, 4 * lg );
	xri[1] = 0.0; /* freq bin at 6400 Hz zeroed */
	ifft3( xri, xnq, (short)lg );
	/*-----------------------------------------------------------*
	* find and quantize gain, multiply xnq[] by gain.           *
	* windowing of xnq[] for TCX overlap.                       *
	*-----------------------------------------------------------*/
	if( pre_echo ) {
		for( i = 0; i < lg; i += 16 ) {
			for( k = 0; k < 16; k++ ) {
				xnq[i + k] *= gain_shap[i / 16];
				xn[i + k] *= gain_shap[i / 16];

				wm[i + k] *= gain_shap[i / 16];
			}
		}
	}
	gain = get_gain( xn, xnq, lg );
	prm[1] = q_gain_tcx( xnq, lg, &gain );
	for( i = 0; i < lg; i++ ) {
		xnq[i] = gain * xnq[i] + gain_pan * wm[i];
	}

	/* adaptive windowing on overlap (beginning and end of frame) */
	for( i = 0; i < ovlp_size; i++ ) {
#ifndef COS_FAC
		xnq[i] *= window[i];
#else
		tmpfloat = cos_fac( i, ovlp_size, lext );
		xnq[i] *= tmpfloat;
#endif
	}
	for( i = 0; i < lext; i++ ) {
#ifndef COS_FAC
		xnq[i + L_frame] *= window[ovlp_size + i];
#else
		tmpfloat = cos_fac( ovlp_size + i, ovlp_size, lext );
		xnq[i + L_frame] *= tmpfloat;
#endif
	}
	for( i = L_frame + lext; i < lg; i++ ) {
		xnq[i] = 0;
	}
	/*-----------------------------------------------------------*
	* TCX overlap and add.  Update memory for next overlap.     *
	*-----------------------------------------------------------*/
	for( i = 0; i < L_OVLP_2k; i++ ) {
		xnq[i] += wovlp[i];
	}
	/* save overlap for next frame */
	for( i = 0; i < lext; i++ ) {
		wovlp[i] = xnq[i + L_frame];
	}
	for( i = lext; i < L_OVLP_2k; i++ ) {
		wovlp[i] = 0.0;
	}
	/*-----------------------------------------------------------*
	* find excitation and synthesis                             *
	*-----------------------------------------------------------*/
	for( i_subfr = 0; i_subfr < L_frame; i_subfr++ ) {
		synth[i_subfr] = xnq[i_subfr];
	}
	return;
}
/*-----------------------------------------------------------------*
* Funtion c_stereo                                             *
* ~~~~~~~~~~~~~~~~~~~~                                            *
*-----------------------------------------------------------------*/
void cod_tcx_stereo( float mono_2k[],
    float                  right_2k[],
    int                    param[],
    int                    brMode,
    int                    mod[],
    int                    fscale,
    Coder_State_Plus *     st )
{
	float sig_buf_mono[TCX_L_FFT_2k];
	float sig_buf_right[TCX_L_FFT_2k];
	float sig_buf_side[TCX_L_FFT_2k];
	float synth_tcx[L_FRAME_2k];
	float synth[L_FRAME_2k];
	int   ovlp_size[4 + 1];
	int   prm_tcx[NPRM_TCX80_D];
	float ovlp[L_OVLP_2k * 5];
	float ovlp_tcx[L_OVLP_2k];
	/* Scalars */
	int   i, k, i2, i1, nbits, *prm;
	float snr, snr1, snr2;
	float tmp;
	/* number of bits per frame (80 ms) */
	nbits = StereoNbits[brMode] - 24 - 4;
	if( StereoNbits[brMode] - 4 > 300 )
		nbits -= 24;

	set_zero( sig_buf_mono, TCX_L_FFT_2k );
	if( fscale == 0 ) {
		mvr2r( mono_2k - TCX_STEREO_DELAY_2k, sig_buf_mono, L_FRAME_2k + L_OVLP_2k );
	}
	else {
		mvr2r( mono_2k - TCX_STEREO_DELAY_2k - L_SUBFR_2k, sig_buf_mono, L_FRAME_2k + L_OVLP_2k );
	}
	set_zero( sig_buf_right, TCX_L_FFT_2k );
	if( fscale == 0 ) {
		mvr2r( right_2k - TCX_STEREO_DELAY_2k, sig_buf_right, L_FRAME_2k + L_OVLP_2k );
	}
	else {
		float *pt = right_2k - TCX_STEREO_DELAY_2k - L_SUBFR_2k;
		mvr2r( right_2k - TCX_STEREO_DELAY_2k - L_SUBFR_2k, sig_buf_right, L_FRAME_2k + L_OVLP_2k );
	}
	set_zero( sig_buf_side, TCX_L_FFT_2k );
	for( i = 0; i < TCX_L_FFT_2k; i++ ) {
		sig_buf_side[i] = sig_buf_mono[i] - sig_buf_right[i];
	}
	/*---------------------------------------------------------------*
	*  Call TCX codec												*
	*---------------------------------------------------------------*/
	ovlp_size[0] = st->mem_stereo_ovlp_size;
	mvr2r( st->mem_stereo_ovlp, ovlp, L_OVLP_2k );
	snr2 = 0.0;
	for( i1 = 0; i1 < 2; i1++ ) {
		snr1 = 0.0;
		for( i2 = 0; i2 < 2; i2++ ) {
			k = ( i1 * 2 ) + i2;
			/* set pointer to parameters */
			prm = param + ( k * NPRM_DIV_TCX_STEREO );
			/*--------------------------------------------------*
			* Call 20MS TCX with pre-echo coder and find segmental SNR *
			*--------------------------------------------------*/
			mvr2r( &ovlp[k * L_OVLP_2k], ovlp_tcx, L_OVLP_2k );
			ctcx_stereo(
			    &sig_buf_side[k * L_DIV_2k],
			    &sig_buf_mono[k * L_DIV_2k],
			    &synth_tcx[k * L_DIV_2k],
			    ovlp_tcx,
			    ovlp_size[k],
			    L_FRAME_2k / 4,
			    nbits / 4 - 2,
			    prm_tcx,
			    1 );
			tmp = segsnr( &sig_buf_side[k * L_DIV_2k], &synth_tcx[k * L_DIV_2k], L_FRAME_2k / 4, L_DIV_2k );
			snr = tmp;
			mod[k] = 0;
			ovlp_size[k + 1] = L_OVLP_2k / 4;
			mvr2r( ovlp_tcx, &ovlp[( k + 1 ) * L_OVLP_2k], L_OVLP_2k );
			mvr2r( &synth_tcx[k * L_DIV_2k], &synth[k * L_DIV_2k], L_FRAME_2k / 4 );
			mvi2i( prm_tcx, prm, NPRM_TCX20_D );
			/*--------------------------------------------------*
			* Call 20MS TCX coder and find segmental SNR       *
			*--------------------------------------------------*/
			mvr2r( &ovlp[k * L_OVLP_2k], ovlp_tcx, L_OVLP_2k );
			ctcx_stereo(
			    &sig_buf_side[k * L_DIV_2k],
			    &sig_buf_mono[k * L_DIV_2k],
			    &synth_tcx[k * L_DIV_2k],
			    ovlp_tcx,
			    ovlp_size[k],
			    L_FRAME_2k / 4,
			    nbits / 4 - 2,
			    prm_tcx,
			    0 );
			tmp = segsnr( &sig_buf_side[k * L_DIV_2k], &synth_tcx[k * L_DIV_2k], L_FRAME_2k / 4, L_DIV_2k );
			if( tmp > snr ) {
				snr = tmp;
				mod[k] = 1;
				ovlp_size[k + 1] = L_OVLP_2k / 4;
				mvr2r( ovlp_tcx, &ovlp[( k + 1 ) * L_OVLP_2k], L_OVLP_2k );
				mvr2r( &synth_tcx[k * L_DIV_2k], &synth[k * L_DIV_2k], L_FRAME_2k / 4 );
				mvi2i( prm_tcx, prm, NPRM_TCX20_D );
			}
			snr1 += 0.5f * snr;
		} /* end of i2 */
		k = ( i1 * 2 );
		/* set pointer to parameters */
		prm = param + ( k * NPRM_DIV_TCX_STEREO );
		/*--------------------------------------------------*
		* Call 40MS TCX coder and find segmental SNR       *
		*--------------------------------------------------*/
		mvr2r( &ovlp[k * L_OVLP_2k], ovlp_tcx, L_OVLP_2k );
		ctcx_stereo(
		    &sig_buf_side[k * L_DIV_2k],
		    &sig_buf_mono[k * L_DIV_2k],
		    &synth_tcx[k * L_DIV_2k],
		    ovlp_tcx,
		    ovlp_size[k],
		    L_FRAME_2k / 2,
		    nbits / 2 - 4,
		    prm_tcx,
		    0 );
		tmp = segsnr( &sig_buf_side[k * L_DIV_2k], &synth_tcx[k * L_DIV_2k], L_FRAME_2k / 2, L_DIV_2k );
		/*--------------------------------------------------------*
		* Save tcx parameters if tcx segmental SNR is better     *
		*--------------------------------------------------------*/
		if( tmp > snr1 ) {
			snr1 = tmp;
			for( i = 0; i < 2; i++ ) {
				mod[k + i] = 2;
			}
			ovlp_size[k + 2] = L_OVLP_2k / 2;
			mvr2r( ovlp_tcx, &ovlp[( k + 2 ) * L_OVLP_2k], L_OVLP_2k );
			mvr2r( &synth_tcx[k * L_DIV_2k], &synth[k * L_DIV_2k], L_FRAME_2k / 2 );
			mvi2i( prm_tcx, prm, NPRM_TCX40_D );
		}
		snr2 += 0.5f * snr1;
	}
	k = 0;
	/* set pointer to parameters */
	prm = param + ( k * NPRM_DIV_TCX_STEREO );
	/*--------------------------------------------------*
	* Call 80MS TCX coder and find segmental SNR       *
	*--------------------------------------------------*/
	mvr2r( &ovlp[k * L_OVLP_2k], ovlp_tcx, L_OVLP_2k );
	ctcx_stereo(
	    sig_buf_side,
	    sig_buf_mono,
	    synth_tcx,
	    ovlp_tcx,
	    ovlp_size[k],
	    L_FRAME_2k,
	    nbits - 8,
	    prm_tcx,
	    0 );
	tmp = segsnr( sig_buf_side, synth_tcx, L_FRAME_2k, L_DIV_2k );
	/*--------------------------------------------------------*
	* Save tcx parameters if tcx segmental SNR is better     *
	*--------------------------------------------------------*/
	if( tmp > snr2 ) {
		snr2 = tmp;
		for( i = 0; i < 4; i++ ) {
			mod[k + i] = 3;
		}
		ovlp_size[k + 4] = L_OVLP_2k;
		mvr2r( ovlp_tcx, &ovlp[( k + 4 ) * L_OVLP_2k], L_OVLP_2k );
		mvr2r( synth_tcx, &synth[k * L_DIV_2k], L_FRAME_2k );
		mvi2i( prm_tcx, prm, NPRM_TCX80_D );
	}
	/*--------------------------------------------------*
	* Update memory.		                            *
	*--------------------------------------------------*/
	st->mem_stereo_ovlp_size = ovlp_size[4];
	mvr2r( &ovlp[4 * L_OVLP_2k], st->mem_stereo_ovlp, L_OVLP_2k );
}
