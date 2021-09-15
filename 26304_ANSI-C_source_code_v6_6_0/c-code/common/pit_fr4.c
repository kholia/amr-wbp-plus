/*-------------------------------------------------------------------*
 * procedure pitch_fr4                                               *
 * ~~~~~~~~~~~~~~~~~~~                                               *
 * Find the closed loop pitch period with 1/4 subsample resolution.  *
 *-------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <math.h>
/* locals functions */
/*-------------------------------------------------------------------*
 * Function  pred_lt4:                                               *
 *           ~~~~~~~~~                                               *
 *-------------------------------------------------------------------*
 * Compute the result of long term prediction with fractionnal       *
 * interpolation of resolution 1/4.                                  *
 *                                                                   *
 * On return exc[0..L_subfr-1] contains the interpolated signal      *
 *   (adaptive codebook excitation)                                  *
 *-------------------------------------------------------------------*/
void pred_lt4(
    float exc[],  /* in/out: excitation buffer */
    int   T0,     /* input : integer pitch lag */
    int   frac,   /* input : fraction of lag   */
    int   L_subfr /* input : subframe size     */
)
{
	int   i, j;
	float s, *x0, *x1, *x2, *c1, *c2;
	x0 = &exc[-T0];
	frac = -frac;
	if( frac < 0 ) {
		frac += PIT_UP_SAMP;
		x0--;
	}
	for( j = 0; j < L_subfr; j++ ) {
		x1 = x0++;
		x2 = x1 + 1;
		c1 = (float *)&inter4_2[frac];
		c2 = (float *)&inter4_2[PIT_UP_SAMP - frac];
		s = 0.0;
		for( i = 0; i < PIT_L_INTERPOL2; i++, c1 += PIT_UP_SAMP, c2 += PIT_UP_SAMP ) {
			s += ( *x1-- ) * ( *c1 ) + ( *x2++ ) * ( *c2 );
		}
		exc[j] = s;
	}
	return;
}
