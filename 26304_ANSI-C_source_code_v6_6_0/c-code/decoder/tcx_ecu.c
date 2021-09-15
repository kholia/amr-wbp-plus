#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*
	phasor computations
*/
void get_phasor( float xri[], float ph[] )
{
	float tmp = 0.0;
	tmp = xri[0] * xri[0] + xri[1] * xri[1];
	if( tmp == 0.0f ) {
		ph[0] = 1;
		ph[1] = 0;
	}
	else {
		tmp = 1.0f / (float)sqrt( tmp );
		ph[0] = xri[0] * tmp;
		ph[1] = xri[1] * tmp;
	}
}
void mult_phasor( float ph1[], float ph2[], float res[] )
{
	float ph_re, ph_im;
	ph_re = ph1[0] * ph2[0] - ph1[1] * ph2[1];
	ph_im = ph1[1] * ph2[0] + ph1[0] * ph2[1];
	res[0] = ph_re;
	res[1] = ph_im;
}
void div_phasor( float ph1[], float ph2[], float res[] )
{
	float ph_re, ph_im;
	ph_re = ph1[0] * ph2[0] + ph1[1] * ph2[1];
	ph_im = ph1[1] * ph2[0] - ph1[0] * ph2[1];
	res[0] = ph_re;
	res[1] = ph_im;
}
/*---------------------------------------------------------------*/
/* Adaptive low frequency de-emphasis in case of packet loss	 */
/*																 */
/*																 */
/*---------------------------------------------------------------*/
void adapt_low_freq_deemph_ecu( float xri[],
    int                               lg,
    Decoder_State_Plus *              st )
{
	int   i, j;
	float buf[L_TCX / 4];
	float max, fac, tmp;
	float pred_max;
	int   curr_mode = 3;
	/*---------------------------------------------------------------*
	* Set correct buffer lengths  								    *
	*---------------------------------------------------------------*/
	if( lg == 1152 ) {
		curr_mode = 3;
	}
	if( lg == 576 ) {
		curr_mode = 2;
	}
	if( lg == 288 ) {
		curr_mode = 1;
	}
	if( ( st->last_mode != curr_mode ) || ( curr_mode <= 2 ) ) {
		pred_max = 0.0f;
	}
	else {
		/*---------------------------------------------------------------*
		* Temporary working buffer									    *
		*---------------------------------------------------------------*/
		mvr2r( st->old_xri, buf, lg / 4 );
		/*---------------------------------------------------------------*
		* Find spectral peak under 1600Hz (lg = 6400Hz)                 *
		* (find maximum of energy of all Re8 subvector under 1600Hz)    *
		*---------------------------------------------------------------*/
		max = 0.01f;
		for( i = 0; i < lg / 4; i += 8 ) {
			tmp = 0.01f;
			for( j = i; j < i + 8; j++ ) {
				tmp += buf[j] * buf[j];
			}
			if( tmp > max ) {
				max = tmp;
			}
		}
		pred_max = (float)sqrt( max ); /* sqrt of energy */
	}
	/*---------------------------------------------------------------*
	* Find spectral peak under 1600Hz (lg = 6400Hz)                 *
	* (find maximum of energy of all Re8 subvector under 1600Hz)    *
	*---------------------------------------------------------------*/
	max = 0.01f;
	for( i = 0; i < lg / 4; i += 8 ) {
		tmp = 0.01f;
		for( j = i; j < i + 8; j++ ) {
			tmp += xri[j] * xri[j];
		}
		if( tmp > max ) {
			max = tmp;
		}
	}
	max = (float)sqrt( max ); /* sqrt of energy */
	/*---------------------------------------------------------------*/
	/* Set the new max												 */
	/*---------------------------------------------------------------*/
	if( max < pred_max ) {
		max = pred_max;
	}
	/*---------------------------------------------------------------*
	* Deemphasis of all subvector below 1600 Hz.                    *
	*---------------------------------------------------------------*/
	fac = 0.1f;
	for( i = 0; i < lg / 4; i += 8 ) {
		tmp = 0.01f;
		for( j = i; j < i + 8; j++ ) {
			tmp += xri[j] * xri[j];
		}
		tmp = (float)sqrt( tmp );
		pessimize();
		tmp = tmp / max;
		if( tmp > fac ) {
			fac = tmp;
		}
		for( j = i; j < i + 8; j++ ) {
			xri[j] *= fac;
		}
	}
}
/*---------------------------------------------------------------*/
/* Spectral reconstruction in case of packet loss				 */
/*																 */
/*																 */
/*---------------------------------------------------------------*/
void reconst_spect( float xri[],
    float                 old_xri[],
    int                   n_pack,
    int                   bfi[],
    int                   lg,
    int                   last_mode,
    float                 buf[] )
{
	float *sp = buf, *old_sp = buf + ( L_TCX / 2 );
	int    lost[L_TCX / 2];
	int    i, l, p;
	int    curr_mode = 3;
	float  energy, old_energy, gain;
	float  gain_sq;
	float  angle, ph_c[2], ph_tmp[2], ph_d[2];
	int    bin, bin_start, bin_end;
	/*---------------------------------------------------------------*
	* Set correct buffer lengths  								    *
	*---------------------------------------------------------------*/
	if( lg == 1152 ) {
		curr_mode = 3;
	}
	if( lg == 576 ) {
		curr_mode = 2;
	}
	if( lg == 288 ) {
		curr_mode = 1;
	}
	/*---------------------------------------------------------------*
	* Can't do anything with ACELP as previous mode  				 *
	*---------------------------------------------------------------*/
	if( ( last_mode != curr_mode ) || ( curr_mode <= 2 ) ) {
		return;
	}
	/*---------------------------------------------------------------*
	* Fill up the lost frequency bins indicator  					 *
	*---------------------------------------------------------------*/
	for( i = 0; i < lg / 2; i++ ) {
		lost[i] = 0;
	}
	for( p = 0; p < n_pack; p++ ) {
		if( bfi[p] ) {
			for( l = p; l < lg / 8; l += n_pack ) {
				/* subvector l is lost */
				for( i = 0; i < 4; i++ ) {
					lost[4 * l + i] = 1;
				}
			}
		}
	}
	/*---------------------------------------------------------------*
	* Compute the old spectrum	       								 *
	*---------------------------------------------------------------*/
	old_sp[0] = old_xri[0] * old_xri[0];
	for( i = 1; i < lg / 2; i++ ) {
		old_sp[i] = old_xri[2 * i] * old_xri[2 * i] + old_xri[2 * i + 1] * old_xri[2 * i + 1];
	}
	/*---------------------------------------------------------------*
	* Compute the new spectrum	       								 *
	*---------------------------------------------------------------*/
	sp[0] = xri[0] * xri[0];
	for( i = 1; i < lg / 2; i++ ) {
		sp[i] = xri[2 * i] * xri[2 * i] + xri[2 * i + 1] * xri[2 * i + 1];
	}
	/*---------------------------------------------------------------*
	* Compute the overall gain difference on non-zero spectral		*
	*---------------------------------------------------------------*/
	energy = 0.0f;
	old_energy = 0.01f;
	for( i = 0; i < lg / 2; i++ ) {
		if( sp[i] > 0.0f ) {
			energy += sp[i];
			old_energy += old_sp[i];
		}
	}
	gain_sq = energy / old_energy;
	pessimize();
	gain = (float)sqrt( gain_sq );
	/* limit the gain */
	if( gain > 1.4142f ) {
		gain = 1.4142f;
		gain_sq = 2.0f;
	}
	/*---------------------------------------------------------------*
	* merge with the new spectrum									*
	*---------------------------------------------------------------*/
	for( i = 0; i < lg; i++ ) {
		if( lost[i / 2] ) {
			xri[i] = gain * old_xri[i];
		}
	}
	/*---------------------------------------------------------------*
	* compensate for the  spectrum									*
	*---------------------------------------------------------------*/
	for( i = 0; i < lg / 2; i++ ) {
		if( lost[i] ) {
			sp[i] = gain_sq * old_sp[i];
		}
	}
	/*---------------------------------------------------------------*
	* apply the group delay conservation							*
	*---------------------------------------------------------------*/
	/* start phase compensation */
	bin = 0;
	do {
		if( lost[bin] ) {
			bin_start = bin;
			while( lost[bin] && bin < lg / 2 ) {
				bin += 4;
			}
			bin_end = bin;
			if( bin_start == 0 ) {
				bin_start = 1;
			}
			/* compute the phase correction factor */
			if( bin_end != lg / 2 ) {
				get_phasor( &xri[2 * ( bin_start - 1 )], ph_tmp );
				get_phasor( &xri[2 * bin_end], ph_c );
				div_phasor( ph_c, ph_tmp, ph_c );
				get_phasor( &old_xri[2 * ( bin_start - 1 )], ph_tmp );
				mult_phasor( ph_c, ph_tmp, ph_c );
				get_phasor( &old_xri[2 * bin_end], ph_tmp );
				div_phasor( ph_c, ph_tmp, ph_c );
				/* get the angle on 4 quad and divide by the length */
				angle = (float)atan2( ph_c[1], ph_c[0] );
				angle /= (float)( bin_end - bin_start + 1 );
				pessimize();
				ph_c[0] = (float)cos( angle );
				ph_c[1] = (float)sin( angle );
			}
			else {
				ph_c[0] = 1;
				ph_c[1] = 0;
			}
			/* loop in the missing bins */
			for( i = bin_start; i < bin_end; i++ ) {
				/* compute old difference */
				get_phasor( &old_xri[2 * ( i - 1 )], ph_tmp );
				get_phasor( &old_xri[2 * i], ph_d );
				div_phasor( ph_d, ph_tmp, ph_d );
				/* compute new phase*/
				get_phasor( &xri[2 * ( i - 1 )], ph_tmp );
				mult_phasor( ph_tmp, ph_d, ph_tmp );
				/* apply phase compensation */
				mult_phasor( ph_tmp, ph_c, ph_tmp );
				xri[2 * i] = (float)sqrt( sp[i] ) * ph_tmp[0];
				xri[2 * i + 1] = (float)sqrt( sp[i] ) * ph_tmp[1];
			}
		}
		bin += 4;
	} while( bin < lg / 2 );
}
