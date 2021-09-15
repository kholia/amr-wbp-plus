/*-------------------------------------------------------------------------*
 * procedure q_gain2_plus                                                  *
 * ~~~~~~~~~~~~~~~~~~~~~~                                                  *
 * Quantization of pitch and codebook gains.                               *
 * The following routines is Q_gains updated for AMR_WB_PLUS.              *
 * MA prediction is removed and MEAN_ENER is now quantized with 2 bits and *
 * transmitted once every ACELP frame to the gains decoder.                *
 * The pitch gain and the code gain are vector quantized and the           *
 * mean-squared weighted error criterion is used in the quantizer search.  *
 *-------------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#define RANGE 64
int q_gain2_plus(     /* (o)  : index of quantizer                      */
    float  code[],    /* (i)  : Innovative code vector                  */
    int    lcode,     /* (i)  : Subframe size                           */
    float *gain_pit,  /* (i/o): Pitch gain / Quantized pitch gain       */
    float *gain_code, /* (i/o): code gain / Quantized codebook gain     */
    float *coeff,     /* (i)  : correlations <y1,y1>, -2<xn,y1>,        */
                      /*                <y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float  mean_ener, /* (i)  : mean_ener defined in open-loop (2 bits) */
    float *g0         /* (o)  : 'correction factor'                     */
)
{
	int          i, indice = 0, min_ind, size;
	float        ener_code, gcode0;
	float        dist, dist_min, g_pitch, g_code;
	const float *t_qua_gain, *p;
	/*-----------------------------------------------------------------*
 * - Find the initial quantization pitch index                     *
 * - Set gains search range                                        *
 *-----------------------------------------------------------------*/
	t_qua_gain = E_ROM_qua_gain7b;
	p = (const float *)( E_ROM_qua_gain7b + RANGE ); /* pt at 1/4th of table */
	min_ind = 0;
	g_pitch = *gain_pit;
	for( i = 0; i < ( NB_QUA_GAIN7B - RANGE ); i++, p += 2 ) {
		if( g_pitch > *p ) {
			min_ind++;
		}
	}
	size = RANGE;
	min_ind = 0;
	size = 128;
	/* innovation energy (without gain) */
	ener_code = 0.01F;
	for( i = 0; i < lcode; i++ ) {
		ener_code += code[i] * code[i];
	}
	ener_code = (float)( 10.0 * log10( ener_code / (float)lcode ) );
	/* predicted codebook gain */
	/* mean energy quantized with 2 bits : 18, 30, 42 or 54 dB */
	gcode0 = mean_ener - ener_code;
	gcode0 = (float)pow( 10.0, gcode0 / 20.0 ); /* predicted gain */
	/* Search for best quantizer */
	dist_min = FLT_MAX;
	p = (const float *)( t_qua_gain + min_ind * 2 );
	for( i = 0; i < size; i++ ) {
		g_pitch = *p++;         /* pitch gain */
		g_code = gcode0 * *p++; /* codebook gain */
		dist = g_pitch * g_pitch * coeff[0]
		    + g_pitch * coeff[1]
		    + g_code * g_code * coeff[2]
		    + g_code * coeff[3]
		    + g_pitch * g_code * coeff[4];
		if( dist < dist_min ) {
			dist_min = dist;
			indice = i;
		}
	}
	indice += min_ind;
	*gain_pit = t_qua_gain[indice * 2];
	*gain_code = t_qua_gain[indice * 2 + 1] * gcode0;
	*g0 = t_qua_gain[indice * 2 + 1];
	return indice;
}
