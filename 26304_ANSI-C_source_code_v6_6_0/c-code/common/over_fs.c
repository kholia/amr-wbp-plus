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

int over_fs(         /* number of sample oversampled       */
    float sig_in[],  /* input:  signal to oversample       */
    float sig_out[], /* output: signal oversampled         */
    int   lg,        /* input:  length of output           */
    int   fac_down,  /* input:  fs*12/fac_down = 44k/48k   */
    float mem[],     /* in/out: mem[2*L_FILT_OVER_FS]      */
    int * frac_mem   /* in/out: interpol fraction memory   */
)
{
	int    i, ncoef, L_frame, L_frame_int;
	float  signal[( 2 * L_FILT_OVER_FS ) + 2 * L_FRAME_MAX];
	float *filter, gain;
	int    fac_up;

	if( fac_down >= 12 )
		return ( lg );

	fac_up = 12;
	filter = (float *)filter_LP12;

	ncoef = L_FILT_OVER_FS;
	gain = 1.0f;
	L_frame = ( ( lg * fac_down ) + ( *frac_mem ) ) / fac_up;
	L_frame_int = lg;

	/* load buffer & update memory */

	for( i = 0; i < ( 2 * L_FILT_OVER_FS ); i++ ) {
		signal[i] = mem[i];
	}
	for( i = 0; i < L_frame; i++ ) {
		signal[i + ( 2 * L_FILT_OVER_FS )] = sig_in[i];
	}

	for( i = 0; i < ( 2 * L_FILT_OVER_FS ); i++ ) {
		mem[i] = signal[i + L_frame];
	}

	/* oversample to 44.1/48 khz */

	interpol_mem( signal + ncoef - 1, sig_out, L_frame_int, filter, ncoef, fac_up, fac_down, gain, frac_mem );

	return ( L_frame );
}

int decim_fs(        /* number of sample decimated         */
    float sig_in[],  /* input:  signal to decimate         */
    int   lg,        /* input:  length of input            */
    float sig_out[], /* output: signal decimated           */
    int   fac_up,    /* input:  44k/48k *fac_up/12 = fs    */
    float mem[],     /* in/out: mem[2*L_FILT_DECIM_FS]     */
    int * frac_mem   /* in/out: interpol fraction memory   */
)
{
	int    i, ncoef, L_frame, L_frame_int;
	float  signal[( 2 * L_FILT_DECIM_FS ) + 2 * L_FRAME_MAX];
	float *filter, gain;
	int    fac_down;

	if( fac_up >= 12 )
		return ( lg );

	fac_down = 12;
	filter = (float *)filter_LP12;

	ncoef = L_FILT_OVER_FS * fac_down / fac_up;
	gain = ( (float)fac_up ) / ( (float)fac_down );
	L_frame = lg;
	L_frame_int = ( ( lg * fac_up ) - ( *frac_mem ) + ( fac_down - 1 ) ) / fac_down;

	/* load buffer & update memory */

	for( i = 0; i < ( 2 * L_FILT_DECIM_FS ); i++ ) {
		signal[i] = mem[i];
	}

	for( i = 0; i < L_frame; i++ ) {
		signal[i + ( 2 * L_FILT_DECIM_FS )] = sig_in[i];
	}

	for( i = 0; i < ( 2 * L_FILT_DECIM_FS ); i++ ) {
		mem[i] = signal[i + L_frame];
	}

	/* decimation from 44.1/48khz to fs_output */

	interpol_mem( signal + ncoef - 1, sig_out, L_frame_int, filter, ncoef, fac_up, fac_down, gain, frac_mem );

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
	int   pos_step_plus_one, fac_up_minus_frac_step;
	float s, *x1, *x2, *c1, *c2;

	pos_step = fac_down / fac_up;
	pos_step_plus_one = pos_step + 1;

	frac_step = fac_down - ( pos_step * fac_up );
	fac_up_minus_frac_step = fac_up - frac_step;

	pos = 0;
	frac = *mem_frac;

	for( i = 0; i < L_frame_int; i++ ) {
		x1 = signal + pos;
		x2 = x1 + 1;
		c1 = &filter[frac];
		c2 = &filter[fac_up - frac];
		s = 0.0;

		for( j = 0; j < nb_coef; j++, c1 += fac_up, c2 += fac_up ) {
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
