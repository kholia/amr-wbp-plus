/*------------------------------------------------------------------*
 * function rnd_ph16                                                *
 * ~~~~~~~~~~~~~~~~~                                                *
 * phase random generator (16 phasef, output lg/2 real+imag)         *
 *------------------------------------------------------------------*/
#include "../include/amr_plus.h"

void rnd_ph16( short *seed, float *xri, int lg )
{
	int            i;
	unsigned short phase;
	for( i = 0; i < lg; i += 2 ) {
		/* random phase from 0 to 15 */
		phase = ( (unsigned short)E_UTIL_random( seed ) ) >> 12;
		xri[i] = sin20[phase + 4];
		xri[i + 1] = sin20[phase];
	}
	return;
}
