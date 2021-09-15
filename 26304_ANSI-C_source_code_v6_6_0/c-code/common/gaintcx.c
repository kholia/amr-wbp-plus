#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
int q_gain_tcx(    /* output: return quantization index */
    float  code[], /* (i)  : quantized vector           */
    int    lcode,  /* (i)  : frame size                 */
    float *gain    /* in/out: quantized gain            */
)
{
	int   i, index;
	float tmp, gcode0, gcode;
	/* energy */
	tmp = 0.01f;
	for( i = 0; i < lcode; i++ ) {
		tmp += code[i] * code[i];
	}
	gcode0 = 4.0f * (float)sqrt( tmp / ( (float)lcode ) );
	/* quantize gain of code */
	gcode = *gain * gcode0;
	if( gcode < FLT_MIN ) {
		gcode = FLT_MIN;
	}
	tmp = 28.0f * (float)log10( gcode ); /* step of 0.714 dB */
	index = (int)floor( tmp + 0.5 );
	if( index < 0 ) {
		index = 0;
	}
	if( index > 127 ) {
		index = 127;
	}
	gcode = (float)pow( 10.0, ( (float)index ) / 28.0 ) / gcode0;
	*gain = gcode;
	return ( index );
}
float d_gain_tcx(    /* output: gain                    */
    int    index,    /* (i)  : index                    */
    float  code[],   /* (i)  : quantized vector         */
    int    lcode,    /* (i)  : frame size               */
    int    bfi,      /* (i)  : 1=gain lost              */
    float *old_rms ) /* (i/o): for frame recovery       */
{
	int   i;
	float tmp, gain, gcode0, rms;
	/* energy */
	tmp = FLT_MIN;
	for( i = 0; i < lcode; i++ ) {
		tmp += code[i] * code[i];
	}
	rms = (float)sqrt( tmp / (float)lcode );
	if( bfi == 0 ) {
		gcode0 = 4.0f * rms;
		gain = (float)pow( 10.0, ( (float)index ) / 28.0 ) / gcode0;
		*old_rms = gain * rms; /* rms of gain*code[] */
	}
	else {
		( *old_rms ) *= .7f; /* decrease energy by -3dB */
		gain = *old_rms / rms;
	}
	return ( gain );
}
