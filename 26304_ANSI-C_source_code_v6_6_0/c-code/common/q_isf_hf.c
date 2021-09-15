/*-------------------------------------------------------------------*
 * functions: q_isf_hf()                                             *
 *                                                                   *
 * Coding/Decoding of HF ISF parameters with AR prediction.          *
 *                                                                   *
 * The ISF vector is quantized using one-stage VQ with 8 elements.   *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#define ISF_GAP_HF 180.0f
#define MU 0.5f
static int sub_VQ( float *x, const float *dico, int dim, int dico_size, float *distance );
void       q_isf_hf(
          float *      isf1,   /* input : ISF in the frequency domain (0..6400) */
          float *      isf_q,  /* output: quantized ISF                         */
          float *      past_q, /* i/o   : past quantized isf (for AR prediction)*/
          int *        indice, /* output: quantization indices (7 words)        */
          const float *mean_isf_hf,
          const float *dico1_isf_hf )
{
	int   i, k, index;
	float min_err, distance, isf2[Q_ISF_ORDER];
	distance = 1.0e30f;
	for( k = 0; k < SIZE_BK1_HF; k++ ) {
		for( i = 0; i < Q_ISF_ORDER; i++ ) {
			isf2[i] = isf1[i] - dico1_isf_hf[( k * Q_ISF_ORDER ) + i] - mean_isf_hf[i] - MU * past_q[i];
		}
		index = sub_VQ( isf2, dico2_isf_hf, Q_ISF_ORDER, SIZE_BK2_HF, &min_err );
		if( min_err < distance ) {
			distance = min_err;
			indice[0] = k;
			indice[1] = index;
		}
	}
	for( i = 0; i < Q_ISF_ORDER; i++ ) {
		isf_q[i] = dico2_isf_hf[( indice[1] * Q_ISF_ORDER ) + i] + dico1_isf_hf[( indice[0] * Q_ISF_ORDER ) + i]
		    + mean_isf_hf[i] + MU * past_q[i];
	}
	E_LPC_isf_reorderPlus( isf_q, ISF_GAP_HF, Q_ISF_ORDER );
	for( i = 0; i < Q_ISF_ORDER; i++ ) {
		past_q[i] = isf_q[i] - mean_isf_hf[i];
	}
	return;
}
/*-------------------------------------------------------------------*
 * routine:   d_isf_hf()                                             *
 *            ~~~~~~~~~                                              *
 * Decoding of ISF parameters                                        *
 *-------------------------------------------------------------------*
 *  Arguments:                                                       *
 *    indice[] : indices of the two selected codebook entries        *
 *    isf[]    : quantized ISFs (in frequency domain)                *
 *-------------------------------------------------------------------*/
#define ALPHA 0.9f
void d_isf_hf(
    int *        indice, /* input:  quantization indices                  */
    float *      isf_q,  /* output: quantized ISFs in the cosine domain   */
    float *      past_q, /* i/o   : past quantized isf (for AR prediction)*/
    int          bfi,    /* input : Bad frame indicator                  */
    const float *mean_isf_hf,
    const float *dico1_isf_hf )
{
	int i;
	if( bfi == 0 ) /* Good frame */
	{
		for( i = 0; i < Q_ISF_ORDER; i++ ) {
			isf_q[i] = dico1_isf_hf[( indice[0] * Q_ISF_ORDER ) + i]
			    + dico2_isf_hf[( indice[1] * Q_ISF_ORDER ) + i]
			    + mean_isf_hf[i] + MU * past_q[i];
		}
	}
	else /* bad frame */
	{
		/* use the past ISFs slightly shifted towards the means */
		for( i = 0; i < Q_ISF_ORDER; i++ ) {
			isf_q[i] = ALPHA * past_q[i] + mean_isf_hf[i];
		}
	}
	E_LPC_isf_reorderPlus( isf_q, ISF_GAP_HF, Q_ISF_ORDER );
	for( i = 0; i < Q_ISF_ORDER; i++ ) {
		past_q[i] = isf_q[i] - mean_isf_hf[i];
	}
	return;
}
/*---------------------------------------------------------------------------*
 * procedure  sub_VQ                                                         *
 *            ~~~~~~                                                         *
 * Quantization of a subvector of size 2 in Split-VQ of ISFs                 *
 *                                                                           *
 * Return:  quantization index                                               *
 *--------------------------------------------------------------------------*/
static int sub_VQ( float *x, const float *dico, int dim, int dico_size, float *distance )
{
	float        dist_min, dist, temp;
	const float *p_dico;
	int          i, j, index;
	dist_min = 1.0e30f;
	p_dico = dico;
	for( i = 0; i < dico_size; i++ ) {
		dist = 0.0;
		for( j = 0; j < dim; j++ ) {
			temp = x[j] - *p_dico++;
			dist += temp * temp;
		}
		if( dist < dist_min ) {
			dist_min = dist;
			index = i;
		}
	}
	*distance = dist_min;
	/* Reading the selected vector */
	p_dico = &dico[index * dim];
	for( j = 0; j < dim; j++ ) {
		x[j] = *p_dico++;
	}
	return index;
}
