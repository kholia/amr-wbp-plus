/*-------------------------------------------------------------------*
 * Function decim_12k8() and oversamp_12k8()                         *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         *
 * decim_12k8    : decimation from fs to 12.8kHz.                    *
 * oversamp_12k8 : oversampling from 12.8kHz to fs.                  *
 *                                                                   *
 * fs = 16/22/24/28.8/32/44/48 kHz.                                  *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define L_FRAME_MAX ( L_FRAME24k )
#define L_FILT_MAX ( L_FILT24k )
#define L_SUBFR16k 80 /* Subframe size at 16kHz           */
static void interpol_add( float *signal, float *signal_int, int L_frame_int, const float *filter, int nb_coef, int fac_up, int fac_down, float gain );
void        decim_12k8(
           float sig_fs[],  /* input:  signal to decimate      */
           int   lg,        /* input:  length of input         */
           float sig12k8[], /* output: decimated signal        */
           float mem[],     /* in/out: memory (2*L_FILT_FS)    */
           int   band )       /* input:  0=0..6.4k, 1=6.4..10.8k */
{
	int          i, ncoef;
	float        signal[( 2 * L_FILT_MAX ) + L_FRAME_MAX];
	const float *filter;
	float        gain;
	int          fac_up, fac_down;
	gain = 1.0;
	switch( lg ) {
	case L_FRAME16kPLUS:
		ncoef = NCOEF_32k;
		fac_up = FAC1_32k;
		fac_down = FAC2_32k;
		if( band == 0 ) {
			filter = filter_32k;
		}
		else {
			if( lg == L_FRAME16kPLUS ) {
				filter = filter_32k_7k;
			}
			else {
				filter = filter_32k_hf;
			}
		}
		break;
	case L_FRAME24k:
		ncoef = NCOEF_48k;
		fac_up = FAC1_48k;
		fac_down = FAC2_48k;
		if( band == 0 ) {
			filter = filter_48k;
		}
		else {
			filter = filter_48k_hf;
		}
		break;
	case L_FRAME8k: /* This mode actually upsamples rather than down */
		ncoef = NCOEF_32k >> 2;
		fac_up = FAC1_32k * 2;
		fac_down = FAC2_32k;
		filter = filter_8k;
		gain = 2.0f;
		break;
	default:
		printf( "wrong frame size in decim_12k8!\n" );
		exit( 0 );
		break;
	}
	if( lg <= L_FRAME24k ) {
		ncoef = ( ncoef + 1 ) >> 1;
		fac_up *= 2;
		gain *= 2.0;
	}
	for( i = 0; i < ( 2 * ncoef ); i++ ) {
		signal[i] = mem[i];
	}
	for( i = 0; i < lg; i++ ) {
		signal[i + ( 2 * ncoef )] = sig_fs[i];
	}
	interpol( signal + ncoef, sig12k8, L_FRAME_PLUS, filter, ncoef, fac_up, fac_down, gain );
	for( i = 0; i < ( 2 * ncoef ); i++ ) {
		mem[i] = signal[i + lg];
	}
	return;
}
void oversamp_12k8(
    float sig12k8[], /* input:  signal to oversampling  */
    float sig_fs[],  /* output: oversampled signal      */
    int   lg,        /* input:  length of output        */
    float mem[],     /* in/out: memory (2*L_FILT)       */
    int   band,      /* input:  0/2/3=0..6.4k, 1=6.4..10.8k */
    int   add        /* input:  1 if result is added to  */
                     /*         output vector, 0 if not */
)
{
	int          i, ncoef, L_frame;
	float        signal[( 2 * L_FILT ) + L_FRAME_PLUS];
	const float *filter;
	float        gain;
	int          fac_up, fac_down;
	gain = 1.0;
	switch( lg ) {
	case L_FRAME16kPLUS:
	case L_SUBFR16k:
		fac_up = FAC2_32k;
		fac_down = FAC1_32k;
		if( band == 0 ) {
			filter = filter_32k;
		}
		else {
			if( lg == L_FRAME16kPLUS ) {
				filter = filter_32k_7k;
			}
			else {
				filter = filter_32k_hf;
			}
		}
		break;
	case L_FRAME24k:
		fac_up = FAC2_48k;
		fac_down = FAC1_48k;
		if( band == 0 ) {
			filter = filter_48k;
		}
		else {
			filter = filter_48k_hf;
		}
		break;
		/* this mode is special: it's actually downsampling, thus needs other filter */
	case L_FRAME8k:
		fac_up = FAC2_32k;
		fac_down = FAC1_32k * 2; /* note: will be multiplied by 2 again below */
		filter = filter_8k;
		gain = 2.0f;
		break;
	default:
		printf( "wrong frame size in oversamp_12k8! (Unsupported rate)\n" );
		exit( 0 );
		break;
	}
	gain *= (float)fac_up / (float)fac_down;
	if( lg <= L_FRAME24k ) {
		fac_down *= 2;
	}
	ncoef = NCOEF_12k8;
	if( lg == L_SUBFR16k )
		L_frame = L_SUBFR;
	else
		L_frame = L_FRAME_PLUS;
	for( i = 0; i < ( 2 * ncoef ); i++ ) {
		signal[i] = mem[i];
	}
	for( i = 0; i < L_frame; i++ ) {
		signal[i + ( 2 * ncoef )] = sig12k8[i];
	}
	if( add == 0 ) {
		interpol( signal + ncoef, sig_fs, lg, filter, ncoef, fac_up, fac_down, gain );
	}
	else {
		interpol_add( signal + ncoef, sig_fs, lg, filter, ncoef, fac_up, fac_down, gain );
	}
	for( i = 0; i < ( 2 * ncoef ); i++ ) {
		mem[i] = signal[i + L_frame];
	}
	return;
}
void interpol( float *signal, float *signal_int, int L_frame_int, const float *filter, int nb_coef, int fac_up, int fac_down, float gain )
{
	int          i, j, frac, frac_step, pos, pos_step;
	int          pos_step_plus_one, fac_up_minus_frac_step;
	float        s, *x1, *x2;
	const float *c1, *c2;
	pos_step = fac_down / fac_up;
	pos_step_plus_one = pos_step + 1;
	frac_step = fac_down - ( pos_step * fac_up );
	fac_up_minus_frac_step = fac_up - frac_step;
	pos = 0;
	frac = 0;
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
		if( frac > fac_up_minus_frac_step ) {
			pos += pos_step_plus_one;
			frac -= fac_up_minus_frac_step;
		}
		else {
			pos += pos_step;
			frac += frac_step;
		}
	}
	return;
}
static void interpol_add( float *signal, float *signal_int, int L_frame_int, const float *filter, int nb_coef, int fac_up, int fac_down, float gain )
{
	int          i, j, frac, frac_step, pos, pos_step;
	int          pos_step_plus_one, fac_up_minus_frac_step;
	float        s, *x1, *x2, tmp;
	const float *c1, *c2;
	pos_step = fac_down / fac_up;
	pos_step_plus_one = pos_step + 1;
	frac_step = fac_down - ( pos_step * fac_up );
	fac_up_minus_frac_step = fac_up - frac_step;
	pos = 0;
	frac = 0;
	for( i = 0; i < L_frame_int; i++ ) {
		x1 = signal + pos;
		x2 = x1 + 1;
		c1 = &filter[frac];
		c2 = &filter[fac_up - frac];
		s = 0.0;
		for( j = 0; j < nb_coef; j++, c1 += fac_up, c2 += fac_up ) {
			s += ( *x1-- ) * ( *c1 ) + ( *x2++ ) * ( *c2 );
		}
		tmp = gain * s; /* this construct must be kept for bit-exactness*/
		signal_int[i] = tmp + signal_int[i];
		if( frac > fac_up_minus_frac_step ) {
			pos += pos_step_plus_one;
			frac -= fac_up_minus_frac_step;
		}
		else {
			pos += pos_step;
			frac += frac_step;
		}
	}
	return;
}
