#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* prototypes*/
void dec_prm_tcx_stereo( int mod[],    /* (i) : frame mode (mode[4], 4 division) */
    int                      param[],  /* (i) : parameters                       */
    short                    serial[], /* (o) : serial bits stream               */
    int                      nbits );

void init_tcx_stereo_decoder( Decoder_State_Plus *st )
{
	st->mem_stereo_ovlp_size = 0;
	set_zero( st->mem_stereo_ovlp, L_OVLP_2k );
	st->mem_balance = 0;
}
float d_balance( /* output: gain                    */
    int index    /* (i)  : index                    */
)
{
	float balance = (float)index / 32.0f - 2.0f;
	return ( balance );
}
void dtcx_stereo( float synth[], /* in/out: synth[-M..lg]           */
    float               mono[],
    float               wovlp[],   /* i/o:    wovlp[0..127]           */
    int                 ovlp_size, /* input:  0, 64 or 128 (0=acelp)  */
    int                 L_frame,   /* input:  frame length            */
    int                 prm[],
    int                 pre_echo,
    int                 bad_frame[],
    Decoder_State_Plus *st )
{
	int    i, k, i_subfr, lg;
	float  tmp, gain;
	int    any_bad_frame;
	float  xri[L_TCX_LB];
	float *xnq = xri;
	float  wm_arr[L_TCX_LB + ECU_WIEN_ORD];
	float  window[32 + 32];
	float *wm = &wm_arr[ECU_WIEN_ORD];
	float  balance;
	float  gain_shap[8];
	float  rmm[ECU_WIEN_ORD + 1];
	float  rms[ECU_WIEN_ORD + 1];

	float *h;
	int    lext = 32;
	int    n_pack = 4;

	/*------ set length of overlap (lext) and length of encoded frame (lg) -----*/
	for( i = 0; i < ECU_WIEN_ORD; i++ ) {
		wm_arr[i] = 0;
	}
	h = st->h;
	switch( L_frame ) {
	case 40:
		lext = 8;
		n_pack = 1;
		break;
	case 80:
		lext = 16;
		n_pack = 2;
		break;
	case 160:
		lext = 32;
		n_pack = 4;
		break;
	};
	lg = L_frame + lext;
	any_bad_frame = 0;
	for( i = 0; i < n_pack; i++ ) {
		any_bad_frame |= bad_frame[i];
	}
	/* built window for overlaps section */
	cos_window( window, ovlp_size, lext );
	/* these are already concealed by the avq demux */
	for( i = 0; i < lg; i++ ) {
		xri[i] = (float)prm[i + 2];
	}
	/* only in mode 0,1 */
	/* windowed mono  */
	for( i = 0; i < lg; i++ ) {
		wm[i] = mono[i];
	}

	/* windowing for TCX overlap and correlation */
	for( i = 0; i < ovlp_size; i++ ) {
		wm[i] *= window[i];
	}
	for( i = 0; i < lext; i++ ) {
		wm[L_frame + i] *= window[ovlp_size + i];
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
		for( i = 0; i < lg / 16; i++ ) {
			gain_shap[i] = (float)sqrt( gain_shap[i] * pow( 10.0f, -tmp ) );

			gain_shap[i] = my_min( 2.0, my_max( 0.5, gain_shap[i] ) );
		}
	}

	/*-----------------------------------------------------------*
		* Compute inverse FFT for obtaining xnq[] without noise.    *
		* Coefficients (xri[]) order are                            *
		*    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
		* Note that last FFT element (re[n/2]) is zeroed.           *
		*-----------------------------------------------------------*/
	adap_low_freq_deemph( xri, 4 * lg );
	xri[1] = 0.0f;
	ifft3( xri, xnq, (short)lg );
	tmp = 1.0f;
	if( bad_frame[0] ) {
		/* if this was predictive we could have a work around */
		balance = (float)( 0.9 * st->mem_balance );
		st->mem_balance = balance;
	}
	else {
		/* balance is safe */
		balance = d_balance( prm[0] );
		st->mem_balance = balance;
	}
	if( pre_echo ) {
		for( i = 0; i < lg; i += 16 ) {
			for( k = 0; k < 16; k++ ) {
				xnq[i + k] *= gain_shap[i / 16];
			}
		}
	}
	gain = d_gain_tcx( prm[1], xnq, lg, any_bad_frame, &( st->side_rms ) );
	if( L_frame != 40 || bad_frame[0] == 0 ) {
		for( i = 0; i < lg; i++ ) {
			xnq[i] = gain * xnq[i] + balance * wm[i];
		}
	}
	else {
		for( i = 0; i < lg; i++ ) {
			xnq[i] = 0;
			for( k = 0; k < ECU_WIEN_ORD + 1; k++ ) {
				xnq[i] += h[k] * wm[i - k];
			}
		}
	}

	if( !any_bad_frame ) {
		/* Create rmm and rms */
		crosscorr( wm, wm, rmm, lg, 0, ECU_WIEN_ORD + 1 );
		crosscorr( xnq, wm, rms, lg, 0, ECU_WIEN_ORD + 1 );
		glev_s( h, rmm, rms, ECU_WIEN_ORD + 1 );
	}
	/*-----------------------------------------------------------*
	* find and quantize gain, multiply xnq[] by gain.           *
	* windowing of xnq[] for TCX overlap.                       *
	*-----------------------------------------------------------*/
	/* adaptive windowing on overlap (beginning and end of frame) */
	for( i = 0; i < ovlp_size; i++ ) {
		xnq[i] *= window[i];
	}
	for( i = 0; i < lext; i++ ) {
		xnq[i + L_frame] *= window[ovlp_size + i];
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
* Funtion dec_tcx_stereo                                           *
* ~~~~~~~~~~~~~~~~~~~~                                            *
*-----------------------------------------------------------------*/
void dec_tcx_stereo( float synth_2k[],
    float                  left_2k[],
    float                  right_2k[],
    int                    param[],
    int                    bad_frame[],
    Decoder_State_Plus *   st )
{
	float synth_side[L_FRAME_2k];
	/* Scalars */
	int  i, k, mod[4];
	int  bfi;
	int *prm;
	bfi = bad_frame[0] || bad_frame[1] || bad_frame[2] || bad_frame[3];
	/* get the mode */
	for( k = 0; k < 4; k++ ) {
		mod[k] = param[k];
	}
	/* tcx decoder loop */
	k = 0;
	while( k < NB_DIV ) {
		prm = param + 4 + ( k * NPRM_DIV_TCX_STEREO );
		if( mod[k] == 1 || mod[k] == 0 ) {
			dtcx_stereo(
			    &synth_side[k * L_DIV_2k],
			    &synth_2k[k * L_DIV_2k],
			    st->mem_stereo_ovlp,
			    st->mem_stereo_ovlp_size,
			    L_FRAME_2k / 4,
			    prm,
			    !mod[k],
			    bad_frame + k,
			    st );
			st->mem_stereo_ovlp_size = L_OVLP_2k / 4;
			k++;
		}
		else if( mod[k] == 2 ) {
			dtcx_stereo(
			    &synth_side[k * L_DIV_2k],
			    &synth_2k[k * L_DIV_2k],
			    st->mem_stereo_ovlp,
			    st->mem_stereo_ovlp_size,
			    L_FRAME_2k / 2,
			    prm,
			    0,
			    bad_frame + k,
			    st );
			st->mem_stereo_ovlp_size = L_OVLP_2k / 2;
			k += 2;
		}
		else if( mod[k] == 3 ) {
			dtcx_stereo(
			    &synth_side[k * L_DIV_2k],
			    &synth_2k[k * L_DIV_2k],
			    st->mem_stereo_ovlp,
			    st->mem_stereo_ovlp_size,
			    L_FRAME_2k,
			    prm,
			    0,
			    bad_frame + k,
			    st );
			st->mem_stereo_ovlp_size = L_OVLP_2k;
			k += 4;
		}
		else {
			printf( "\n error, unknown mode" );
			exit( -1 );
		}
	}
	k = 0;
	while( k < NB_DIV ) {

		if( mod[k] == 0 || mod[k] == 1 ) {
			for( i = k * L_FRAME_2k / 4; i < k * L_FRAME_2k / 4 + L_FRAME_2k / 4; i++ ) {
				left_2k[i] = synth_2k[i] + synth_side[i];
				right_2k[i] = synth_2k[i] - synth_side[i];
			}
			k++;
		}
		else if( mod[k] == 2 ) {
			for( i = k * L_FRAME_2k / 4; i < k * L_FRAME_2k / 4 + L_FRAME_2k / 2; i++ ) {
				left_2k[i] = synth_2k[i] + synth_side[i];
				right_2k[i] = synth_2k[i] - synth_side[i];
			}
			k += 2;
		}
		else {
			for( i = 0; i < L_FRAME_2k; i++ ) {
				left_2k[i] = synth_2k[i] + synth_side[i];
				right_2k[i] = synth_2k[i] - synth_side[i];
			}
			k += 4;
		}
	}
	return;
}
