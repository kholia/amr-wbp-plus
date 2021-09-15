/*-------------------------------------------------------------------*
 * functions: qpisf_2s() and dpisf_2s()                              *
 *                                                                   *
 * Coding/Decoding of ISF parameters  with prediction.               *
 *                                                                   *
 * The ISF vector is quantized using two-stage VQ with split-by-2    *
 * in 1st stage and split-by-5 in the second stage.                  *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#define MU ( 1.0f / 3.0f ) /* prediction factor */
#define N_SURV_MAX 4       /* 4 survivors max   */
/*------------------------------------------------------------------*
 * routine:   qpisf_2s()                                            *
 *            ~~~~~~~~~                                             *
 * Quantization of ISF parameters with prediction                   *
 *                                                                  *
 * The ISF vector is quantized using two-stage VQ with split-by-2   *
 * in 1st stage and split-by-5 in the second stage.                 *
 *------------------------------------------------------------------*/
static int  sub_VQ( float *x, const float *dico, int dim, int dico_size, float *err );
static void VQ_stage1( float *x, const float *dico, int dim, int dico_size, int *index, int surv );
void        qpisf_2s_46b(
           float *isf1,      /* input : ISF in the frequency domain (0..6400) */
           float *isf_q,     /* output: quantized ISF                         */
           float *past_isfq, /* i/0   : past ISF quantizer                    */
           int *  indice,    /* output: quantization indices (7 words)        */
           int    nb_surv    /* input : number of survivor (1, 2, 3 or 4)     */
       )
{
	int   i, k, tmp_ind[5];
	int   surv1[N_SURV_MAX]; /* indices of survivors from 1st stage */
	float temp, min_err, distance, isf[ORDER];
	float isf_stage2[ORDER];
	for( i = 0; i < ORDER; i++ ) {
		isf[i] = isf1[i] - mean_isf[i] - MU * past_isfq[i];
	}
	VQ_stage1( &isf[0], dico1_isf, 9, SIZE_BK1, surv1, nb_surv );
	distance = 1.0e30f;
	if( nb_surv > N_SURV_MAX ) {
		nb_surv = N_SURV_MAX;
	}
	for( k = 0; k < nb_surv; k++ ) {
		for( i = 0; i < 9; i++ ) {
			isf_stage2[i] = isf[i] - dico1_isf[i + surv1[k] * 9];
		}
		tmp_ind[0] = sub_VQ( &isf_stage2[0], dico21_isf, 3, SIZE_BK21, &min_err );
		temp = min_err;
		tmp_ind[1] = sub_VQ( &isf_stage2[3], dico22_isf, 3, SIZE_BK22, &min_err );
		temp += min_err;
		tmp_ind[2] = sub_VQ( &isf_stage2[6], dico23_isf, 3, SIZE_BK23, &min_err );
		temp += min_err;
		if( temp < distance ) {
			distance = temp;
			indice[0] = surv1[k];
			for( i = 0; i < 3; i++ ) {
				indice[i + 2] = tmp_ind[i];
			}
		}
	}
	VQ_stage1( &isf[9], dico2_isf, 7, SIZE_BK2, surv1, nb_surv );
	distance = 1.0e30f;
	for( k = 0; k < nb_surv; k++ ) {
		for( i = 0; i < 7; i++ ) {
			isf_stage2[i] = isf[9 + i] - dico2_isf[i + surv1[k] * 7];
		}
		tmp_ind[0] = sub_VQ( &isf_stage2[0], dico24_isf, 3, SIZE_BK24, &min_err );
		temp = min_err;
		tmp_ind[1] = sub_VQ( &isf_stage2[3], dico25_isf, 4, SIZE_BK25, &min_err );
		temp += min_err;
		if( temp < distance ) {
			distance = temp;
			indice[1] = surv1[k];
			for( i = 0; i < 2; i++ ) {
				indice[i + 5] = tmp_ind[i];
			}
		}
	}
	/* decoding the ISFs */
	dpisf_2s_46b( indice, isf_q, past_isfq, isf_q, isf_q, 0, 0, 0 );
	return;
}
/*-------------------------------------------------------------------*
 * routine:   dpisf_2s()                                              *
 *            ~~~~~~~~~                                              *
 * Decoding of ISF parameters                                        *
 *-------------------------------------------------------------------*
 *  Arguments:                                                       *
 *    indice[] : indices of the two selected codebook entries        *
 *    isf[]    : quantized ISFs (in frequency domain)                *
 *-------------------------------------------------------------------*/
#define ALPHA 0.9f
void dpisf_2s_46b(
    int *  indice,     /* input:  quantization indices                  */
    float *isf_q,      /* output: quantized ISFs in the cosine domain   */
    float *past_isfq,  /* i/0   : past ISF quantizer                    */
    float *isfold,     /* input : past quantized ISF                    */
    float *isf_buf,    /* input : isf buffer                            */
    int    bfi,        /* input : Bad frame indicator                   */
    int    bfi_2nd_st, /* input : 2nd stage bfi mask (bin: 011111)      */
    int    enc_dec )
{
	int   i, j;
	float tmp;
	float ref_isf[ORDER];
	if( bfi == 0 ) /* Good frame */
	{
		for( i = 0; i < 9; i++ ) {
			isf_q[i] = dico1_isf[indice[0] * 9 + i];
		}
		for( i = 0; i < 7; i++ ) {
			isf_q[i + 9] = dico2_isf[indice[1] * 7 + i];
		}
		if( ( bfi_2nd_st & 1 ) == 0 ) {
			for( i = 0; i < 3; i++ ) {
				isf_q[i] += dico21_isf[indice[2] * 3 + i];
			}
		}
		if( ( bfi_2nd_st & 2 ) == 0 ) {
			for( i = 0; i < 3; i++ ) {
				isf_q[i + 3] += dico22_isf[indice[3] * 3 + i];
			}
		}
		if( ( bfi_2nd_st & 4 ) == 0 ) {
			for( i = 0; i < 3; i++ ) {
				isf_q[i + 6] += dico23_isf[indice[4] * 3 + i];
			}
		}
		if( ( bfi_2nd_st & 8 ) == 0 ) {
			for( i = 0; i < 3; i++ ) {
				isf_q[i + 9] += dico24_isf[indice[5] * 3 + i];
			}
		}
		if( ( bfi_2nd_st & 16 ) == 0 ) {
			for( i = 0; i < 4; i++ ) {
				isf_q[i + 12] += dico25_isf[indice[6] * 4 + i];
			}
		}
		for( i = 0; i < ORDER; i++ ) {
			tmp = isf_q[i];
			isf_q[i] = tmp + MU * past_isfq[i] + mean_isf[i];
			past_isfq[i] = tmp;
		}
		if( enc_dec ) {
			for( i = 0; i < ORDER; i++ ) {
				for( j = ( L_MEANBUF - 1 ); j > 0; j-- ) {
					isf_buf[j * ORDER + i] = isf_buf[( j - 1 ) * ORDER + i];
				}
				isf_buf[i] = isf_q[i];
			}
		}
	}
	else /* bad frame */
	{
		/* use the past ISFs slightly shifted towards their mean */
		for( i = 0; i < ORDER; i++ ) {
			ref_isf[i] = mean_isf[i] / ( L_MEANBUF + 1 );
			for( j = 0; j < L_MEANBUF; j++ ) {
				ref_isf[i] += isf_buf[j * ORDER + i] / ( L_MEANBUF + 1 );
			}
		}
		for( i = 0; i < ORDER; i++ ) {
			isf_q[i] = ALPHA * isfold[i] + ( 1.0f - ALPHA ) * ref_isf[i];
		}
		/* estimate past quantized residual to be used in next frame */
		for( i = 0; i < ORDER; i++ ) {
			tmp = ref_isf[i] + past_isfq[i] * MU; /* predicted ISF */
			past_isfq[i] = isf_q[i] - tmp;
			past_isfq[i] *= 0.5;
		}
	}
	E_LPC_isf_reorderPlus( isf_q, ISF_GAP, ORDER );
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
	index = 0;
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
static void VQ_stage1( float *x, const float *dico, int dim, int dico_size, int *index, int surv )
{
	float        dist_min[N_SURV_MAX];
	float        dist, temp;
	const float *p_dico;
	int          i, j, k, l;
	for( i = 0; i < surv; i++ ) {
		dist_min[i] = 1.0e30f;
	}
	for( i = 0; i < surv; i++ ) {
		index[i] = i;
	}
	p_dico = dico;
	for( i = 0; i < dico_size; i++ ) {
		dist = 0.0;
		for( j = 0; j < dim; j++ ) {
			temp = x[j] - *p_dico++;
			dist += temp * temp;
		}
		for( k = 0; k < surv; k++ ) {
			if( dist < dist_min[k] ) {
				for( l = surv - 1; l > k; l-- ) {
					dist_min[l] = dist_min[l - 1];
					index[l] = index[l - 1];
				}
				dist_min[k] = dist;
				index[k] = i;
				break;
			}
		}
	}
	return;
}
