#include "../include/amr_plus.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define L_FRAME_MAX ( L_FRAME48k )

/* local function */
static void interpol_mem(
    float *signal,
    float *signal_int,
    int    L_frame_int,
    float *filter,
    int    nb_coef,
    int    fac_up,
    int    fac_down,
    float  gain,
    int *  mem_frac );

int decim_split_12k8(   /* number of sample decimated         */
    float sig_fs[],     /* input:  signal to decimate         */
    int   lg_input,     /* input:  2*L_FRAME44k if 44kHz      */
    float sig12k8_lf[], /* output: LF decimated signal        */
    float sig12k8_hf[], /* output: HF decimated signal        */
    int   lg,           /* input:  length of LF and HF        */
    int   fac_fs,       /* input:  at 48kHz, scale fac = fac_fs/FSCALE_DENOM */
    int   fac_up,       /* input:  Upsampling factor */
    int   fac_down,     /* input:  Downsampling factor */
    int   L_frame,      /* input:  Working frame len */
    float mem[],        /* in/out: mem[L_MEM_DECIM_SPLIT]     */
    int * frac_mem      /* in/out: interpol fraction memory   */
)
{
	int    i, j, ncoef, L_frame_int;
	float  signal[( 2 * ( L_FILT_SPLIT + L_FILT_DECIM ) ) + 2 * L_FRAME_MAX], *sig;
	float *filter, gain, *x1, *x2, s;

#ifdef FILTER_44kHz
	if( lg_input == ( 2 * L_FRAME44k ) ) {
		filter = (float *)filter_LP165;
	}
#endif
#ifdef FILTER_48kHz
	if( lg_input == ( 2 * L_FRAME48k ) ) {
		filter = (float *)filter_LP180;
	}
#endif

	ncoef = ( L_FILT_OVER + 1 ) * fac_down / fac_up;
	gain = ( (float)fac_up ) / ( (float)fac_down );
	/* frames length */
	L_frame_int = 2 * lg; /* 25k6 rate */
	/* load buffer & update memory */
	for( i = 0; i < ( 2 * ( L_FILT_SPLIT + L_FILT_DECIM ) ); i++ ) {
		signal[i] = mem[i];
	}
	for( i = 0; i < L_frame; i++ ) {
		signal[i + ( 2 * ( L_FILT_SPLIT + L_FILT_DECIM ) )] = sig_fs[i];
	}

	sig = signal + ( 2 * L_FILT_SPLIT );
	for( i = 0; i < ( 2 * L_FILT_DECIM ); i++ ) {
		mem[i + ( 2 * L_FILT_SPLIT )] = sig[i + L_frame];
	}
	/* decimation to 25.6kHz */
	interpol_mem( sig + ncoef, sig, L_frame_int, filter, ncoef, fac_up, fac_down, gain, frac_mem );
	/* update memory */
	for( i = 0; i < ( 2 * L_FILT_SPLIT ); i++ ) {
		mem[i] = signal[i + L_frame_int];
	}
	/* band split (25.6kHz --> 2 x 12.8kHz) */
	sig = signal + L_FILT_SPLIT;
	for( i = 0; i < lg; i++ ) {
		s = sig[i * 2];
		x1 = &sig[( i * 2 ) - 1];
		x2 = &sig[( i * 2 ) + 1];
		for( j = 0; j < L_FILT_JOIN; j++, x1 -= 2, x2 += 2 ) {
			s += ( *x1 + *x2 ) * inter2_coef[j];
		}
		s *= 0.5f;
		sig12k8_lf[i] = s;
		sig12k8_hf[i] = sig[i * 2] - s;
	}
	return ( L_frame );
}
int join_over_12k8(     /* number of sample oversampled       */
    float sig12k8_lf[], /* input:  LF signal (fs=12k8)        */
    float sig12k8_hf[], /* input:  HF signal (fs=12k8)        */
    int   lg,           /* input:  length of LF and HF        */
    float sig_fs[],     /* output: oversampled signal         */
    int   lg_output,    /* input:  L_FRAME44k if 44kHz        */
    int   fac_fs,       /* input:  at 48kHz, scale fac = fac_fs/FSCALE_DENOM */
    float mem[],        /* in/out: mem[L_MEM_JOIN_OVER]       */
    int * frac_mem      /* in/out: interpol fraction memory   */
)
{
	int    i, j, ncoef, L_frame, L_frame_int;
	float  signal[( 2 * L_FILT_JOIN ) + L_FRAME_PLUS];
	float  signal2[( 2 * L_FILT_OVER ) + ( 2 * L_FRAME_PLUS )], *sig;
	float *filter, gain, *x1, *x2, s;
	int    fac_up, fac_down;

#ifdef FILTER_44kHz
	if( lg_output == L_FRAME44k ) {
		fac_up = 3 * 441;
		fac_down = fac_fs * 8;
		filter = (float *)filter_LP165;
	}
#endif
#ifdef FILTER_48kHz
	if( lg_output == L_FRAME48k ) {
		fac_up = 180 * 8;
		fac_down = fac_fs * 8;
		filter = (float *)filter_LP180;
	}
#endif

	ncoef = L_FILT_OVER;
	gain = 1.0f;
	/* frames length */
	L_frame = 2 * lg; /* 25k6 rate */
	L_frame_int = ( ( L_frame * fac_up ) - ( *frac_mem ) + ( fac_down - 1 ) ) / fac_down;
	/* band join (2 x 12.8kHz --> 25.6kHz) */
	/* load buffer (LF-HF) & update memory */
	for( i = 0; i < ( 2 * L_FILT_JOIN ); i++ ) {
		signal[i] = mem[i + ( 2 * L_FILT_OVER )];
	}
	for( i = 0; i < lg; i++ ) {
		signal[i + ( 2 * L_FILT_JOIN )] = sig12k8_lf[i] - sig12k8_hf[i];
	}
	for( i = 0; i < ( 2 * L_FILT_JOIN ); i++ ) {
		mem[i + ( 2 * L_FILT_OVER )] = signal[i + lg];
	}
	/* odd samples = LF-HF interpolated */
	sig = signal2 + ( 2 * L_FILT_OVER );
	for( i = 0; i < lg; i++ ) {
		x1 = &signal[i + L_FILT_JOIN];
		x2 = &signal[i + 1 + L_FILT_JOIN];
		s = 0.0;
		for( j = 0; j < L_FILT_JOIN; j++, x1--, x2++ ) {
			s += ( *x1 + *x2 ) * inter2_coef[j];
		}
		sig[( i * 2 ) + 1] = s;
	}
	/* synchronise LF+HF with interpolated LF-HF */
	for( i = 0; i < ( 2 * L_FILT_JOIN ); i++ ) {
		signal[i] = mem[i + ( 2 * L_FILT_OVER ) + ( 2 * L_FILT_JOIN )];
	}
	for( i = 0; i < lg; i++ ) {
		signal[i + ( 2 * L_FILT_JOIN )] = sig12k8_lf[i] + sig12k8_hf[i];
	}
	for( i = 0; i < ( 2 * L_FILT_JOIN ); i++ ) {
		mem[i + ( 2 * L_FILT_OVER ) + ( 2 * L_FILT_JOIN )] = signal[i + lg];
	}
	/* even samples = LF+HF */
	for( i = 0; i < lg; i++ ) {
		sig[i * 2] = signal[i + L_FILT_JOIN];
	}
	/* load buffer & update memory (for oversampling) */
	for( i = 0; i < ( 2 * L_FILT_OVER ); i++ ) {
		signal2[i] = mem[i];
	}
	for( i = 0; i < ( 2 * L_FILT_OVER ); i++ ) {
		mem[i] = signal2[i + L_frame];
	}
	/* oversample from 25.6kHz to fs */
	interpol_mem( signal2 + ncoef, sig_fs, L_frame_int, filter, ncoef, fac_up, fac_down, gain, frac_mem );
	return ( L_frame_int );
}
static void interpol_mem(
    float *signal,
    float *signal_int,
    int    L_frame_int,
    float *filter,
    int    nb_coef,

    int   fac_up,
    int   fac_down,
    float gain,
    int * mem_frac )
{
	int   i, j, frac, frac_step, pos, pos_step;
	int   frac_d, fac_up_d;
	int   pos_step_plus_one, fac_up_minus_frac_step;
	float s, *x1, *x2, *c1, *c2;

	pos_step = fac_down / fac_up;
	pos_step_plus_one = pos_step + 1;

	frac_step = fac_down - ( pos_step * fac_up );
	fac_up_minus_frac_step = fac_up - frac_step;

	fac_up_d = fac_up >> 3;
	pos = 0;
	frac = *mem_frac;

	for( i = 0; i < L_frame_int; i++ ) {
		x1 = signal + pos;
		x2 = x1 + 1;
		frac_d = frac >> 3;
		c1 = &filter[frac_d];
		c2 = &filter[fac_up_d - frac_d];

		s = 0.0;

		for( j = 0; j < nb_coef; j++, c1 += fac_up_d, c2 += fac_up_d ) {
			s += ( *x1-- ) * ( *c1 ) + ( *x2++ ) * ( *c2 );
		}
		signal_int[i] = gain * s;

		if( frac >= fac_up_minus_frac_step ) {
			pos += pos_step_plus_one;
			frac -= fac_up_minus_frac_step;
		}
		else {
			pos += pos_step;
			frac += frac_step;
		}
	}
	*mem_frac = frac;
	return;
}
