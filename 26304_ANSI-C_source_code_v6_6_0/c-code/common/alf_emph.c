#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
/*---------------------------------------------------------------*
 * Adaptive low frequencies emphasis (for SV below 1600 hz).     *
 *                                                               *
 * To ensure that low frequencies is well quantized (otherwise   *
 * audible overlaps problem occurs), each subvector (RE8) is     *
 * is emphased relatively to the spectral peak under 1600Hz.     *
 *                                                               *
 * The subvector gain (dB) is:                                   *
 *    0.5 * min( spectral peak (dB) / local peak (dB), 40dB)     *
 * where                                                         *
 *    spectral peak = peak between 0 Hz and 1600 Hz.             *
 *    local peak = peak between 0Hz to the subvector.            *
 *                                                               *
 * The gain is never over 20dB at the beginning and decrease to  *
 * 0 dB somewhere between 0 and 1600 Hz.                         *
 *---------------------------------------------------------------*/
void adap_low_freq_emph( float xri[], int lg )
{
	int   i, j;
	float max, fac, tmp;
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
	                          /*---------------------------------------------------------------*
  * Emphasis of all subvector below 1600 Hz.                      *
  *---------------------------------------------------------------*/
	fac = 10.0f;
	for( i = 0; i < lg / 4; i += 8 ) {
		tmp = 0.01f;
		for( j = i; j < i + 8; j++ ) {
			tmp += xri[j] * xri[j];
		}
		tmp = (float)sqrt( tmp );
		tmp = (float)sqrt( max / tmp );
		if( tmp < fac ) {
			fac = tmp;
		}
		for( j = i; j < i + 8; j++ ) {
			xri[j] *= fac;
		}
	}
	return;
}
/*---------------------------------------------------------------*
 * Adaptive low frequencies weak deemphasis (for SV below 1600hz)*
 *                                                               *
 * To ensure that low frequencies is well quantized (otherwise   *
 * audible overlaps problem occurs), each subvector (RE8) was    *
 * emphased relatively to the spectral peak under 1600Hz.        *
 *                                                               *
 * This routine is the inverse operation of low_freq_emph().     *
 *---------------------------------------------------------------*/
void adap_low_freq_deemph( float xri[], int lg )
{
	int   i, j;
	float max, fac, tmp;
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
	return;
}
