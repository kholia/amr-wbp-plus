#include <math.h>
#define PI2 6.283185307
#include "../include/amr_plus.h"
/*----------------------------------------------------------*
 * Procedure  cos_window:                                   *
 *            ~~~~~~~~~~~                                   *
 *    To find the cos window of length n1+n2                *
 *                                                          *
 * fh[i] = cos(-pi/2 ... pi/2)                              *
 *----------------------------------------------------------*/
void cos_window( float *fh, int n1, int n2 )
{
	double cc, cte;
	int    i;
	cte = 0.25 * PI2 / (float)n1;
	cc = 0.5 * cte - 0.25 * PI2;
	for( i = 0; i < n1; i++ ) {
		*fh++ = (float)cos( cc );
		cc += cte;
	}
	cte = 0.25 * PI2 / (float)n2;
	cc = 0.5 * cte;
	for( i = 0; i < n2; i++ ) {
		*fh++ = (float)cos( cc );
		cc += cte;
	}
	return;
}
