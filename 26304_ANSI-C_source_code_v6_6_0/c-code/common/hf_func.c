#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
/*---------------------------------------------------------------------*
 * routine match_gain_6k4()                                            *
 * ~~~~~~~~~~~~~~~~~~~~~~~~                                            *
 * Find the gain to match amplitude of both filters at 6.4kHz.         *
 * Gain is returned in log domain.                                     *
 *---------------------------------------------------------------------*/
float match_gain_6k4( float *AqLF, float *AqHF )
{
	float tmp, gain;
	float buf[M + L_SUBFR];
	float code[L_SUBFR];
	int   i;
	/* generate 1/(1+0.9 z^{-1}) over one subframe (HF signal close to a sinusoid at 6.4 kHz) */
	set_zero( buf, M );
	tmp = 1.0;
	for( i = 0; i < L_SUBFR; i++ ) {
		buf[i + M] = tmp;
		tmp *= -0.9f;
	}
	/* apply A(z) from lower-band and 1/A(z) from upper-band */
	E_UTIL_residu( AqLF, buf + M, code, L_SUBFR );
	E_UTIL_synthesisPlus( AqHF, MHF, code, code, L_SUBFR, buf, 0 );
	/* compute energy and correction scale factor */
	tmp = 0.0;
	for( i = 0; i < L_SUBFR; i++ ) {
		tmp += code[i] * code[i];
	}
	gain = 10.0f * (float)log10( 1.0f / tmp );
	return ( gain );
}
/*---------------------------------------------------------------------*
 * routine int_gain()                                                  *
 * ~~~~~~~~~~~~~~~~~~                                                  *
 * Find the interpolated gain parameters in every subframes.           *
 *---------------------------------------------------------------------*/
void int_gain( float old_gain, float new_gain, float *gain, int nb_subfr )
{
	float inc, fnew, fold;
	int   k;
	inc = 1.0f / (float)nb_subfr;
	fnew = 0.0f;
	for( k = 0; k < nb_subfr; k++ ) {
		fold = 1.0f - fnew;
		*gain = ( old_gain * fold ) + ( new_gain * fnew );
		fnew += inc;
		gain++;
	}
	return;
}
/*---------------------------------------------------------------------*
 * routine soft_exc_hf()                                               *
 * ~~~~~~~~~~~~~~~~~~~~~                                               *
 * reduce buzziness from excitation.                                   *
 *---------------------------------------------------------------------*/
void soft_exc_hf( float *exc_hf, float *mem )
{
	float tmp, lp_amp;
	int   i;
	lp_amp = *mem;
	for( i = 0; i < L_SUBFR; i++ ) {
		tmp = (float)fabs( exc_hf[i] );
		lp_amp = 0.98f * lp_amp + 0.02f * tmp;
		tmp = tmp - 2.0f * lp_amp;
		if( tmp < 0.0 ) {
			tmp = 0.0;
		}
		lp_amp += 0.5f * tmp;
		if( exc_hf[i] >= 0.0 ) {
			exc_hf[i] -= tmp;
		}
		else {
			exc_hf[i] += tmp;
		}
	}
	*mem = lp_amp;
	return;
}
/*---------------------------------------------------------------------*
 * routine soft_exc_hf_new()                                           *
 * ~~~~~~~~~~~~~~~~~~~~~                                               *
 * reduce buzziness from excitation.                                   *
 *---------------------------------------------------------------------*/
void soft_exc_hf_new( float *exc_hf, float *mem, int l_frame )
{
	float tmp, lp_amp;
	int   i;
	lp_amp = *mem;
	for( i = 0; i < l_frame; i++ ) {
		tmp = (float)fabs( exc_hf[i] );
		lp_amp = 0.98f * lp_amp + 0.02f * tmp;
		tmp = tmp - 2.0f * lp_amp;
		if( tmp < 0.0 ) {
			tmp = 0.0;
		}
		lp_amp += 0.5f * tmp;
		if( exc_hf[i] >= 0.0 ) {
			exc_hf[i] -= tmp;
		}
		else {
			exc_hf[i] += tmp;
		}
	}
	*mem = lp_amp;
	return;
}
/*---------------------------------------------------------------------*
 * routine smooth_ener_hf()                                            *
 * ~~~~~~~~~~~~~~~~~~~~~~~~                                            *
 * smooth energy evolution of HF synthesis subframe.                   *
 *---------------------------------------------------------------------*/
void smooth_ener_hf( float *HF, float *threshold )
{
	float tmp, ener;
	int   i;
	/* compute energy over subframe */
	ener = 0.0001f;
	for( i = 0; i < L_SUBFR; i++ ) {
		ener += HF[i] * HF[i];
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
	tmp = (float)sqrt( tmp / ener );
	for( i = 0; i < L_SUBFR; i++ ) {
		HF[i] *= tmp;
	}
	return;
}
