#include "../include/amr_plus.h"
#include "../include/s_util.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void init_dec_hi_stereo( Decoder_State_Plus *st )
{
	int i;
	cos_window( st->w_window, 0, L_SUBFR );
	for( i = 0; i < L_SUBFR; i++ ) {
		st->w_window[i] = st->w_window[i] * st->w_window[i];
	}
	set_zero( st->old_wh, HI_FILT_ORDER );
	set_zero( st->old_wh2, HI_FILT_ORDER );
	set_zero( st->old_exc_mono, HI_FILT_ORDER );
	set_zero( st->old_AqLF, 5 * ( M + 1 ) );
	for( i = 0; i < 5; i++ ) {
		st->old_AqLF[i * ( M + 1 )] = 1.0f;
	}
	set_zero( st->left.mem_synth_hi, M );
	set_zero( st->right.mem_synth_hi, M );
	set_zero( st->old_wh_q, HI_FILT_ORDER );
	set_zero( st->old_gm_gain, 2 );
	for( i = 0; i < 4; i++ ) {
		st->old_gain_left[i] = 1.0f;
		st->old_gain_right[i] = 1.0f;
	}
}
static void dec_filt( int **prm,
    float                   wh[],
    float                   old_wh[],
    int                     bfi,
    const PMSVQ *           filt_hi_pmsvq )
{
	pmsvq_inv( wh, prm, old_wh, bfi, filt_hi_pmsvq );
}
static void dec_gain( int **prm,
    float *                 gain_left,
    float *                 gain_right,
    float *                 old_gain,
    int                     bfi,
    const PMSVQ *           gain_hi_pmsvq )
{
	float tmp[2];
	pmsvq_inv( tmp, prm, old_gain, bfi, gain_hi_pmsvq );
	*gain_left = tmp[0];
	*gain_right = tmp[1];
}
/* this supposes time alignement between encoder and decoder */
void dec_hi_stereo( float synth_hi_t0[],
    float                 right_hi[],
    float                 left_hi[],
    float                 AqLF[],
    int                   param[],
    int                   bad_frame[],
    int                   fscale,
    Decoder_State_Plus *  st )
{
	float *exc_buf = right_hi;
	float *exc_mono = exc_buf + HI_FILT_ORDER;
	float *synth_hi_d;
	float *p_Aq, *p_h;
	int    i, k;
	int    i_subfr;
	int    start;
	float  gain_left[NB_DIV + 2];
	float  gain_right[NB_DIV + 2];
	float  new_gain_left[NB_DIV];
	float  new_gain_right[NB_DIV];
	float *left_exc, *right_exc;
	float  side_buf[L_DIV + L_SUBFR];
	float  buf[L_SUBFR];
	float  wh[NB_DIV * HI_FILT_ORDER];
	int *  prm;
	/*---------------------------------------------------------------------*
	* Mismatch between received parameters and signal buffer              *
	*                                                                     *
	*                                                                     *
	*                                                                     *
	*	fscale=0														  *
	*                                                                     *
	*                  |    20ms   |    20ms   |     20ms  |   20ms    |  *
	*        -|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|  *
	*        <--------  Total speech buffer ---------------->			  *
	*                                                                     *
	*	fscale!=0                                                         *
	*                                                                     *
	*                  |    20ms   |    20ms   |     20ms  |   20ms    |  *
	*     -|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|  *
	*     <--------  Total speech buffer ---------------->                *
	* g(-2)    g(-1)       g(0)         g(1)        g(2)        g(3)      *
	* f(-2)    f(-1)       f(0)         f(1)        f(2)        f(3)      *
	*---------------------------------------------------------------------*/
	left_exc = left_hi;
	right_exc = right_hi;
	/*-----------------------------------------------------------*
	* Decode gain and filters				    		         *
	*									                         *
	*------------------------------------------------------------*/
	for( i = 0; i < NB_DIV; i++ ) {
		prm = param + i * NPRM_STEREO_HI_X;
		dec_filt( &prm, &wh[i * HI_FILT_ORDER], st->old_wh_q, bad_frame[i], st->filt_hi_pmsvq );
		dec_gain( &prm, &new_gain_left[i], &new_gain_right[i], st->old_gm_gain, bad_frame[i], st->gain_hi_pmsvq );
		new_gain_left[i] = (float)pow( 10.0f, new_gain_left[i] / 20.0f );
		new_gain_right[i] = (float)pow( 10.0f, new_gain_right[i] / 20.0f );
	}
	/*-----------------------------------------------------------*
	* Set correct alignement of parameters	    		         *
	*------------------------------------------------------------*/
	mvr2r( st->old_gain_left, gain_left, 2 );
	mvr2r( st->old_gain_right, gain_right, 2 );
	mvr2r( new_gain_left, gain_left + 2, 4 );
	mvr2r( new_gain_right, gain_right + 2, 4 );
	mvr2r( new_gain_left + 2, st->old_gain_left, 2 );
	mvr2r( new_gain_right + 2, st->old_gain_right, 2 );
	/*-----------------------------------------------------------*
	* Set correct starting buffer					.           *
	*-----------------------------------------------------------*/
	synth_hi_d = &synth_hi_t0[-L_BSP - D_BPF];
	/*-----------------------------------------------------------*
	* Compute left right pseudo-excitation	    		         *
	*Use the mono LPC filter for whitening                       *
	*------------------------------------------------------------*/
	p_Aq = st->old_AqLF;
	residu( p_Aq, M, synth_hi_d, exc_mono, L_SUBFR / 2 );
	p_Aq += ( M + 1 );
	for( i_subfr = L_SUBFR / 2; i_subfr < 3 * L_SUBFR + L_SUBFR / 2; i_subfr += L_SUBFR ) {
		residu( p_Aq, M, &synth_hi_d[i_subfr], &exc_mono[i_subfr], L_SUBFR );
		p_Aq += ( M + 1 );
	}
	if( fscale != 0 ) {
		residu( p_Aq, M, &synth_hi_d[i_subfr], &exc_mono[i_subfr], L_SUBFR );
	}
	p_Aq = AqLF;
	if( fscale == 0 ) {
		start = 3 * L_SUBFR + L_SUBFR / 2;
	}
	else {
		start = 4 * L_SUBFR + L_SUBFR / 2;
	}
	for( i_subfr = start; i_subfr < L_FRAME_PLUS - L_SUBFR / 2; i_subfr += L_SUBFR ) {
		residu( p_Aq, M, &synth_hi_d[i_subfr], &exc_mono[i_subfr], L_SUBFR );
		p_Aq += ( M + 1 );
	}
	residu( p_Aq, M, &synth_hi_d[L_FRAME_PLUS - L_SUBFR / 2], &exc_mono[L_FRAME_PLUS - L_SUBFR / 2], L_SUBFR / 2 );
	/*-----------------------------------------------------------*
	* prepare excitation buffer for filtering	   		         *
	*------------------------------------------------------------*/
	mvr2r( st->old_exc_mono, exc_buf, HI_FILT_ORDER );
	mvr2r( exc_buf + L_FRAME_PLUS, st->old_exc_mono, HI_FILT_ORDER );
	/*-----------------------------------------------------------*
	* synthesise left right pseudo excitation with interpolation *
	*------------------------------------------------------------*/
	if( fscale == 0 ) {
		fir_filt( st->old_wh, HI_FILT_ORDER, exc_mono, side_buf, L_DIV - L_SUBFR / 2 );
		fir_filt( st->old_wh2, HI_FILT_ORDER, exc_mono, buf, L_SUBFR / 2 );
		for( i = 0; i < L_SUBFR / 2; i++ ) {
			left_exc[i] = st->w_window[i + L_SUBFR / 2] * gain_left[0] * ( buf[i] + exc_mono[i] ) + gain_left[1] * ( 1.0f - st->w_window[i + L_SUBFR / 2] ) * ( exc_mono[i] + side_buf[i] );
			right_exc[i] = st->w_window[i + L_SUBFR / 2] * gain_right[0] * ( exc_mono[i] - buf[i] ) + gain_right[1] * ( 1.0f - st->w_window[i + L_SUBFR / 2] ) * ( exc_mono[i] - side_buf[i] );
		}
		for( i = L_SUBFR / 2; i < L_DIV - L_SUBFR / 2; i++ ) {
			left_exc[i] = gain_left[1] * ( exc_mono[i] + side_buf[i] );
			right_exc[i] = gain_right[1] * ( exc_mono[i] - side_buf[i] );
		}
		fir_filt( st->old_wh, HI_FILT_ORDER, &exc_mono[L_DIV - L_SUBFR / 2], buf, L_SUBFR );
		p_h = wh;
		k = 1;
		for( i_subfr = L_DIV - L_SUBFR / 2; i_subfr < L_FRAME_PLUS - L_SUBFR / 2; i_subfr += L_DIV ) {
			fir_filt( p_h, HI_FILT_ORDER, &exc_mono[i_subfr], side_buf, L_DIV );
			for( i = 0; i < L_SUBFR; i++ ) {
				left_exc[i_subfr + i] = st->w_window[i] * gain_left[k] * ( exc_mono[i_subfr + i] + buf[i] ) + gain_left[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[i_subfr + i] + side_buf[i] );
				right_exc[i_subfr + i] = st->w_window[i] * gain_right[k] * ( exc_mono[i_subfr + i] - buf[i] ) + gain_right[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[i_subfr + i] - side_buf[i] );
			}
			for( i = L_SUBFR; i < L_DIV; i++ ) {
				left_exc[i_subfr + i] = gain_left[k + 1] * ( exc_mono[i_subfr + i] + side_buf[i] );
				right_exc[i_subfr + i] = gain_right[k + 1] * ( exc_mono[i_subfr + i] - side_buf[i] );
			}
			if( i_subfr + L_DIV + L_SUBFR < L_FRAME_PLUS - L_SUBFR / 2 ) {
				fir_filt( p_h, HI_FILT_ORDER, &exc_mono[i_subfr + L_DIV], buf, L_SUBFR );
			}
			else {
				fir_filt( p_h, HI_FILT_ORDER, &exc_mono[i_subfr + L_DIV], buf, L_SUBFR / 2 );
			}
			p_h += HI_FILT_ORDER;
			k++;
		}
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[L_FRAME_PLUS - L_SUBFR / 2], side_buf, L_SUBFR / 2 );
		for( i = 0; i < L_SUBFR / 2; i++ ) {
			left_exc[L_FRAME_PLUS - L_SUBFR / 2 + i] = st->w_window[i] * gain_left[k] * ( exc_mono[L_FRAME_PLUS - L_SUBFR / 2 + i] + buf[i] ) + gain_left[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_FRAME_PLUS - L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[L_FRAME_PLUS - L_SUBFR / 2 + i] = st->w_window[i] * gain_right[k] * ( exc_mono[i_subfr + i] - buf[i] ) + gain_right[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_FRAME_PLUS - L_SUBFR / 2 + i] - side_buf[i] );
		}
	}
	else {
		fir_filt( st->old_wh2, HI_FILT_ORDER, exc_mono, buf, L_SUBFR / 2 );
		for( i = 0; i < L_SUBFR / 2; i++ ) {
			left_exc[i] = gain_left[0] * ( exc_mono[i] + buf[i] );
			right_exc[i] = gain_right[0] * ( exc_mono[i] - buf[i] );
		}
		fir_filt( st->old_wh2, HI_FILT_ORDER, exc_mono + L_SUBFR / 2, buf, L_SUBFR );
		fir_filt( st->old_wh, HI_FILT_ORDER, exc_mono + L_SUBFR / 2, side_buf, L_DIV );
		for( i = 0; i < L_SUBFR; i++ ) {
			left_exc[L_SUBFR / 2 + i] = st->w_window[i] * gain_left[0] * ( exc_mono[L_SUBFR / 2 + i] + buf[i] ) + gain_left[1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[L_SUBFR / 2 + i] = st->w_window[i] * gain_right[0] * ( exc_mono[L_SUBFR / 2 + i] - buf[i] ) + gain_right[1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_SUBFR / 2 + i] - side_buf[i] );
		}
		for( i = L_SUBFR + L_SUBFR / 2; i < L_DIV + L_SUBFR / 2; i++ ) {
			left_exc[i] = gain_left[1] * ( exc_mono[i] + side_buf[i - L_SUBFR / 2] );
			right_exc[i] = gain_right[1] * ( exc_mono[i] - side_buf[i - L_SUBFR / 2] );
		}
		fir_filt( st->old_wh, HI_FILT_ORDER, &exc_mono[L_DIV + L_SUBFR / 2], buf, L_SUBFR );
		p_h = wh;
		k = 1;
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[L_DIV + L_SUBFR / 2], side_buf, L_DIV );
		for( i = 0; i < L_SUBFR; i++ ) {
			left_exc[L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_left[k] * ( exc_mono[L_DIV + L_SUBFR / 2 + i] + buf[i] ) + gain_left[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_right[k] * ( exc_mono[L_DIV + L_SUBFR / 2 + i] - buf[i] ) + gain_right[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		for( i = L_SUBFR; i < L_DIV; i++ ) {
			left_exc[L_DIV + L_SUBFR / 2 + i] = gain_left[k + 1] * ( exc_mono[L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[L_DIV + L_SUBFR / 2 + i] = gain_right[k + 1] * ( exc_mono[L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[( 2 * L_DIV ) + ( L_SUBFR / 2 )], buf, L_SUBFR );
		p_h += HI_FILT_ORDER;
		k++;
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[( 2 * L_DIV ) + L_SUBFR / 2], side_buf, L_DIV );
		for( i = 0; i < L_SUBFR; i++ ) {
			left_exc[2 * L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_left[k] * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] + buf[i] ) + gain_left[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[2 * L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_right[k] * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] - buf[i] ) + gain_right[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		for( i = L_SUBFR; i < L_DIV; i++ ) {
			left_exc[2 * L_DIV + L_SUBFR / 2 + i] = gain_left[k + 1] * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[2 * L_DIV + L_SUBFR / 2 + i] = gain_right[k + 1] * ( exc_mono[2 * L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[( 3 * L_DIV ) + ( L_SUBFR / 2 )], buf, L_SUBFR );
		p_h += HI_FILT_ORDER;
		k++;
		fir_filt( p_h, HI_FILT_ORDER, &exc_mono[( 3 * L_DIV ) + L_SUBFR / 2], side_buf, L_DIV - L_SUBFR / 2 );
		for( i = 0; i < L_SUBFR; i++ ) {
			left_exc[3 * L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_left[k] * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] + buf[i] ) + gain_left[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[3 * L_DIV + L_SUBFR / 2 + i] = st->w_window[i] * gain_right[k] * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] - buf[i] ) + gain_right[k + 1] * ( 1.0f - st->w_window[i] ) * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		for( i = L_SUBFR; i < L_DIV - L_SUBFR / 2; i++ ) {
			left_exc[3 * L_DIV + L_SUBFR / 2 + i] = gain_left[k + 1] * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] + side_buf[i] );
			right_exc[3 * L_DIV + L_SUBFR / 2 + i] = gain_right[k + 1] * ( exc_mono[3 * L_DIV + L_SUBFR / 2 + i] - side_buf[i] );
		}
		p_h += HI_FILT_ORDER;
	}
	/*-----------------------------------------------------------*
	* Save filters in memory for next frame						 *
	*------------------------------------------------------------*/
	mvr2r( p_h, st->old_wh, HI_FILT_ORDER );
	mvr2r( p_h - HI_FILT_ORDER, st->old_wh2, HI_FILT_ORDER );
	/*-----------------------------------------------------------*
	*															 *
	* Synthesise left right mid band signals					 *
	*															 *
	*------------------------------------------------------------*/
	p_Aq = st->old_AqLF;
	syn_filt( p_Aq, M, &left_exc[0], left_hi, L_SUBFR / 2, st->left.mem_synth_hi, 1 );
	syn_filt( p_Aq, M, &right_exc[0], right_hi, L_SUBFR / 2, st->right.mem_synth_hi, 1 );
	p_Aq += ( M + 1 );
	for( k = 0; k < 3; k++ ) {
		syn_filt( p_Aq, M, &left_exc[k * L_SUBFR + L_SUBFR / 2], &left_hi[k * L_SUBFR + L_SUBFR / 2], L_SUBFR, st->left.mem_synth_hi, 1 );
		syn_filt( p_Aq, M, &right_exc[k * L_SUBFR + L_SUBFR / 2], &right_hi[k * L_SUBFR + L_SUBFR / 2], L_SUBFR, st->right.mem_synth_hi, 1 );
		p_Aq += ( M + 1 );
	}
	if( fscale != 0 ) {
		syn_filt( p_Aq, M, &left_exc[k * L_SUBFR + L_SUBFR / 2], &left_hi[k * L_SUBFR + L_SUBFR / 2], L_SUBFR, st->left.mem_synth_hi, 1 );
		syn_filt( p_Aq, M, &right_exc[k * L_SUBFR + L_SUBFR / 2], &right_hi[k * L_SUBFR + L_SUBFR / 2], L_SUBFR, st->right.mem_synth_hi, 1 );
	}
	if( fscale == 0 ) {
		start = 3 * L_SUBFR + L_SUBFR / 2;
	}
	else {
		start = 4 * L_SUBFR + L_SUBFR / 2;
	}
	p_Aq = AqLF;
	for( i_subfr = start; i_subfr < L_FRAME_PLUS - L_SUBFR / 2; i_subfr += L_SUBFR ) {
		syn_filt( p_Aq, M, &left_exc[i_subfr], &left_hi[i_subfr], L_SUBFR, st->left.mem_synth_hi, 1 );
		syn_filt( p_Aq, M, &right_exc[i_subfr], &right_hi[i_subfr], L_SUBFR, st->right.mem_synth_hi, 1 );
		p_Aq += ( M + 1 );
	}
	syn_filt( p_Aq, M, &left_exc[L_FRAME_PLUS - L_SUBFR / 2], &left_hi[L_FRAME_PLUS - L_SUBFR / 2], L_SUBFR / 2, st->left.mem_synth_hi, 1 );
	syn_filt( p_Aq, M, &right_exc[L_FRAME_PLUS - L_SUBFR / 2], &right_hi[L_FRAME_PLUS - L_SUBFR / 2], L_SUBFR / 2, st->right.mem_synth_hi, 1 );
	if( fscale == 0 ) {
		mvr2r( p_Aq, st->old_AqLF, 4 * ( M + 1 ) );
	}
	else {
		mvr2r( p_Aq, st->old_AqLF, 5 * ( M + 1 ) );
	}
}
