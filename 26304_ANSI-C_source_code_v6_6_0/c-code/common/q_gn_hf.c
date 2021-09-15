/*-------------------------------------------------------------------*
 * functions: q_gain_hf()                                            *
 *                                                                   *
 * Coding/Decoding of HF gain parameters with AR prediction.         *
 *                                                                   *
 * The gain vector is quantized using one-stage VQ with 4 elements.  *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
void q_gain_hf(
    float *gain,   /* input : gain of 4 subfr */
    float *gain_q, /* output: quantized gains */
    int *  indice  /* output: indices         */
)
{
	int   i, k;
	float min_err, distance, gain2[Q_GN_ORDER];
	distance = 1.0e30f;
	for( k = 0; k < SIZE_BK_HF; k++ ) {
		for( i = 0; i < Q_GN_ORDER; i++ ) {
			gain2[i] = gain[i] - dico_gain_hf[( k * Q_GN_ORDER ) + i] - MEAN_GAIN_HF;
		}
		min_err = 0.0;
		for( i = 0; i < Q_GN_ORDER; i++ ) {
			min_err += gain2[i] * gain2[i];
		}
		if( min_err < distance ) {
			distance = min_err;
			indice[0] = k;
		}
	}
	for( i = 0; i < Q_GN_ORDER; i++ ) {
		gain_q[i] = dico_gain_hf[( indice[0] * Q_GN_ORDER ) + i] + MEAN_GAIN_HF;
	}
	return;
}
/*-------------------------------------------------------------------*
 * routine:   d_gain_hf()                                            *
 *            ~~~~~~~~~                                              *
 * Decoding of gain parameters                                       *
 *-------------------------------------------------------------------*
 *  Arguments:                                                       *
 *    indice[] : indices of the two selected codebook entries        *
 *    gain[]   : quantized gains                                     *
 *-------------------------------------------------------------------*/
#define ALPHA 0.9f
void d_gain_hf(
    int    indice, /* input:  quantization indices */
    float *gain_q, /* output: quantized gains      */
    float *past_q, /* i/o   : past quantized gain (1 word) */
    int    bfi     /* input : Bad frame indicator  */
)
{
	int   i;
	float tmp;
	if( bfi == 0 ) /* Good frame */
	{
		for( i = 0; i < Q_GN_ORDER; i++ ) {
			gain_q[i] = dico_gain_hf[( indice * Q_GN_ORDER ) + i] + MEAN_GAIN_HF;
		}
	}
	else /* bad frame */
	{
		/* use the past gains slightly shifted towards the means */
		*past_q = ( ALPHA * ( *past_q + 20.0f ) ) - 20.0f;
		for( i = 0; i < Q_GN_ORDER; i++ ) {
			gain_q[i] = *past_q + MEAN_GAIN_HF;
		}
	}
	tmp = 0.0;
	for( i = 0; i < Q_GN_ORDER; i++ ) {
		tmp += gain_q[i];
	}
	*past_q = 0.25f * tmp - MEAN_GAIN_HF;
	return;
}
