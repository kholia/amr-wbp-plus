/*-------------------------------------------------------------------*
 * procedure d_gain2_plus                                            *
 * ~~~~~~~~~~~~~~~~~~~~~~                                            *
 * Decoding of pitch and codebook gains  (see q_gain2_plus.c)        *
 *-------------------------------------------------------------------*
 * input arguments:                                                  *
 *                                                                   *
 *   indice     :Quantization index                                  *
 *   code[]     :Innovative code vector                              *
 *   lcode      :Subframe size                                       *
 *   bfi        :Bad frame indicator                                 *
 *                                                                   *
 * output arguments:                                                 *
 *                                                                   *
 *   gain_pit   :Quantized pitch gain                                *
 *   gain_code  :Quantized codeebook gain                            *
 *                                                                   *
 * Global variables defining quantizer (in qua_gns.h)                *
 *                                                                   *
 *   t_qua_gain[]    :Table of gain quantizers                       *
 *   nb_qua_gain     :Nombre de quantization levels                  *
 *                                                                   *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
float d_gain2_plus(     /* (o)  : 'correction factor' */
    int    index,       /* (i)  : index of quantizer                      */
    float  code[],      /* (i)  : Innovative code vector                  */
    int    lcode,       /* (i)  : Subframe size                           */
    float *gain_pit,    /* (o)  : Quantized pitch gain                    */
    float *gain_code,   /* (o)  : Quantized codebook gain                 */
    int    bfi,         /* (i)  : Bad frame indicato                      */
    float  mean_ener,   /* (i)  : mean_ener defined in open-loop (2 bits) */
    float *past_gpit,   /* (i)  : past gain of pitch                      */
    float *past_gcode ) /* (i/o): past gain of code                       */
{
	int          i;
	float        ener_code, gcode0;
	const float *t_qua_gain;
	float        ener_inov, gcode_inov;
	t_qua_gain = E_ROM_qua_gain7b;
	ener_inov = 0.01f;
	for( i = 0; i < lcode; i++ ) {
		ener_inov += code[i] * code[i];
	}
	gcode_inov = (float)( 1.0 / sqrt( ener_inov / (float)lcode ) );
	/*----------------- Test erasure ---------------*/
	if( bfi != 0 ) {
		if( *past_gpit > 0.95f ) {
			*past_gpit = 0.95f;
		}
		if( *past_gpit < 0.5f ) {
			*past_gpit = 0.5f;
		}
		*gain_pit = *past_gpit;
		*past_gpit *= 0.95f;
		*past_gcode *= ( 1.4f - *past_gpit );
		*gain_code = *past_gcode * gcode_inov;
		return *past_gcode;
	}
	/*-------------- Decode gains ---------------*/
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
	*gain_pit = t_qua_gain[index * 2];
	*gain_code = t_qua_gain[index * 2 + 1] * gcode0;
	/* update bad frame handler */
	*past_gpit = *gain_pit;
	*past_gcode = *gain_code / gcode_inov;
	return t_qua_gain[index * 2 + 1];
}
