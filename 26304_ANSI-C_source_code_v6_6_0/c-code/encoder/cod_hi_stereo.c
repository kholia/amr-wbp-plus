#include "../include/amr_plus.h"
#include "../include/s_util.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void init_cod_hi_stereo( Coder_State_Plus *st )
{
	int i;
	set_zero( st->old_exc_mono, HI_FILT_ORDER );
	st->filt_energy_threshold = 0.0f;
	set_zero( st->old_wh, HI_FILT_ORDER );
	set_zero( st->old_wh_q, HI_FILT_ORDER );
	set_zero( st->old_gm_gain, 2 );
	cos_window( st->w_window, 0, L_SUBFR );
	for( i = 0; i < L_SUBFR; i++ ) {
		st->w_window[i] = st->w_window[i] * st->w_window[i];
	}
}
/* gain  quantizer */
static void quant_gain( float gain_left,  /* i/o */
    float                     gain_right, /* i/o */
    float                     old_gain[],
    int **                    prm,
    const PMSVQ *             gain_hi_pmsvq )
{
	float tmp[2];
	tmp[0] = gain_left;
	tmp[1] = gain_right;
	pmsvq( tmp, prm, tmp, old_gain, gain_hi_pmsvq );
}
/* filter  quantizer */
static void quant_filt( float h[],     /* i/o */
    float                     old_h[], /* i/o */
    int **                    prm,
    const PMSVQ *             filt_hi_pmsvq )
{
	pmsvq( h, prm, h, old_h, filt_hi_pmsvq );
}
/* filter smoother */
static void smooth_ener_filter( float *filter,
    float *                            threshold )
{
	float tmp, ener, old_ener;
	int   i;
	/* compute energy over subframe */
	ener = 0.0001f;
	for( i = 0; i < HI_FILT_ORDER; i++ ) {
		ener += filter[i] * filter[i];
	}
	/* in any case limit the filter energy */
	old_ener = ener;
	if( ener > 16.0f ) {
		ener = 16.0f;
	}
	/* if energy<threshold, add 1.5 dB and saturate to threshold
     else substract 1.5 dB and saturate to threshold */
	tmp = ener;
	if( tmp < *threshold ) {
		tmp = tmp * 1.414f;
		if( tmp > *threshold ) {
			tmp = *threshold;
		}
	}
	else {
		tmp = tmp / 1.414f;
		if( tmp < *threshold ) {
			tmp = *threshold;
		}
	}
	/* set the threshold for next subframer to the current modified energy */
	*threshold = tmp;
	/* apply correction scale factor to HF signal */
	tmp = (float)sqrt( tmp / old_ener );
	for( i = 0; i < HI_FILT_ORDER; i++ ) {
		filter[i] *= tmp;
	}
	return;
}
void cod_hi_stereo( float speech_hi[],
    float                 right_hi[],
    float                 AqLF[],
    int                   param[],
    Coder_State_Plus *    st )
{
	float *exc_mono = speech_hi - M;
	float *exc_side = right_hi - M;
	int    i_subfr, i, k, t, j;
	float *p_Aq;
	/* covariance matrix */
	float r[HI_FILT_ORDER][HI_FILT_ORDER];
	float c[HI_FILT_ORDER];
	/* estimated LMS filters */
	float  wh[NB_DIV * HI_FILT_ORDER];
	float *p_h;
	/* signal big subframe pointers */
	float *x, *y;
	float  buf[L_DIV + L_SUBFR];
	float  gain_left[NB_SUBFR];
	float  gain_right[NB_SUBFR];
	/* energies */
	float energy_right;
	float energy_right_q;
	float energy_left;
	float energy_mono;
	float energy_left_q;
	float corr_left_right;
	float gain_fact;
	int * prm;
	mvr2r( st->old_exc_mono, exc_mono - HI_FILT_ORDER, HI_FILT_ORDER );
	/* compute the residual of the hi mono and right */
	p_Aq = AqLF;
	for( i_subfr = 0; i_subfr < L_FRAME_PLUS; i_subfr += L_SUBFR ) {
		residu( p_Aq, M, &speech_hi[i_subfr], &exc_mono[i_subfr], L_SUBFR );
		residu( p_Aq, M, &right_hi[i_subfr], &exc_side[i_subfr], L_SUBFR );
		p_Aq += ( M + 1 );
	}
	residu( p_Aq, M, &speech_hi[i_subfr], &exc_mono[i_subfr], L_SUBFR );
	residu( p_Aq, M, &right_hi[i_subfr], &exc_side[i_subfr], L_SUBFR );
	/* compute real side signal */
	for( i = 0; i < L_FRAME_PLUS + L_SUBFR; i++ ) {
		exc_side[i] = exc_mono[i] - exc_side[i];
	}
	/* save fir state for next frame */
	mvr2r( exc_mono + L_FRAME_PLUS - HI_FILT_ORDER, st->old_exc_mono, HI_FILT_ORDER );
	/* compute the wiener filters, raw on each frame with covariance method*/
	p_h = wh;
	for( i = 0; i < NB_DIV; i++ ) {
		/* set the pointer to parameters */
		prm = param + i * NPRM_STEREO_HI_X;
		/* set signal pointers */
		x = exc_mono + i * L_DIV;
		y = exc_side + i * L_DIV;
		/* compute cross-correlation terms */
		for( k = 0; k < HI_FILT_ORDER; k++ ) {
			c[k] = 0.0f;
			for( t = 0; t < L_DIV; t++ ) {
				c[k] += y[t] * x[t - k];
			}
			for( t = L_DIV; t < L_DIV + L_SUBFR; t++ ) {
				c[k] += y[t] * x[t - k];
			}
		}
		/* compute correlation matrix */
		for( k = 0; k < HI_FILT_ORDER; k++ ) {
			for( j = k; j < HI_FILT_ORDER; j++ ) {
				r[k][j] = 0.0f;
				for( t = 0; t < L_DIV; t++ ) {
					r[k][j] += x[t - k] * x[t - j];
				}
				for( t = L_DIV; t < L_DIV + L_SUBFR; t++ ) {
					r[k][j] += x[t - k] * x[t - j];
				}
			}
		}
		/* compute a solution to the linear system */
		if( cholsolc( r, c, p_h, HI_FILT_ORDER ) ) {
			/* cholesky failed use panning */
			for( k = 1; k < HI_FILT_ORDER; k++ ) {
				p_h[k] = 0.0f;
			}
			p_h[0] = c[0] / ( r[0][0] + 1.0f );
		}
		/* wiener filter energy smoothing  */
		smooth_ener_filter( p_h, &st->filt_energy_threshold );
		/* quantize the filters*/
		quant_filt( p_h, st->old_wh_q, &prm, st->filt_hi_pmsvq );
		/* local synthesis */
		fir_filt( p_h, HI_FILT_ORDER, x, buf, L_DIV + L_SUBFR );
		/* compute the gain matching figures in the excitation domain */
		energy_left = 0.001f;
		energy_right = 0.001f;
		energy_left_q = 0.001f;
		energy_right_q = 0.001f;
		energy_mono = 0.001f;
		corr_left_right = 0.00f;
		for( t = 0; t < L_DIV + L_SUBFR; t++ ) {
			/* mono + side */
			energy_left += ( x[t] + y[t] ) * ( x[t] + y[t] );
			energy_left_q += ( x[t] + buf[t] ) * ( x[t] + buf[t] );
			/* mono */
			energy_mono += x[t] * x[t];
			/* mono - side */
			energy_right += ( x[t] - y[t] ) * ( x[t] - y[t] );
			energy_right_q += ( x[t] - buf[t] ) * ( x[t] - buf[t] );
			corr_left_right += ( x[t] - y[t] ) * ( x[t] + y[t] );
		}
		/* compute the gain matching */
		gain_right[i] = 10.0f * (float)log10( energy_right / energy_right_q + 0.001f );
		gain_left[i] = 10.0f * (float)log10( energy_left / energy_left_q + 0.001f );
		gain_fact = ( 4.0f * energy_mono ) / ( energy_left + energy_right );
		if( gain_fact < 1.0f ) {
			/* we have signal cancelation, take no risks */
			if( gain_right[i] > 0.0f ) {
				gain_right[i] = my_max( gain_right[i] + 10.0f * (float)log10( gain_fact ), 0.0f );
			}
			if( gain_left[i] > 0.0f ) {
				gain_left[i] = my_max( gain_left[i] + 10.0f * (float)log10( gain_fact ), 0.0f );
			}
		}
		/* quantize the gains */
		quant_gain( gain_left[i], gain_right[i], st->old_gm_gain, &prm, st->gain_hi_pmsvq );
		/* next frame */
		p_h += HI_FILT_ORDER;
	}
	/* save last filter*/
	mvr2r( &wh[( NB_DIV - 1 ) * HI_FILT_ORDER], st->old_wh, HI_FILT_ORDER );
}
