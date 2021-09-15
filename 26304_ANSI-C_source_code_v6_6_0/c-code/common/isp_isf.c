/*-----------------------------------------------------------------------*
 *   Functions : isp2isf and isf2isp                                     *
 *                                                                       *
 *      isp2isf   Transformation isp to isf                              *
 *      isf2isp   Transformation isf to isp                              *
 *                                                                       *
 *   ISP are immitance spectral pair in cosine domain (-1 to 1).         *
 *   ISF are immitance spectral pair in frequency domain (0 to 6400).    *
 *-----------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#define SCALE1 ( 6400.0 / PI )
#define SCALE2 ( PI / 6400.0 )
#ifndef PI
#define PI 3.141592654
#endif
void isf2isp(
    float isf[], /* input : isf[m] normalized (range: 0<=val<=6400)  */
    float isp[], /* output: isp[m] (range: -1<=val<1)                */
    int   m      /* input : LPC order                                */
)
{
	int i;
	/*  convert ISFs to the cosine domain */
	for( i = 0; i < m - 1; i++ ) {
		isp[i] = (float)cos( isf[i] * SCALE2 );
	}
	isp[m - 1] = (float)cos( isf[m - 1] * SCALE2 * 2.0 );
	return;
}
