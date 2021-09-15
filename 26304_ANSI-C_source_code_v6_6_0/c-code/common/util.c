#include "../include/amr_plus.h"
#include <math.h>
/*-----------------------------------------------------------*
 * procedure set_zero                                        *
 * ~~~~~~~~~~~~~~~~~~                                        *
 * Set a vector x[] of dimension n to zero.                  *
 *-----------------------------------------------------------*/
void set_zero( float *x, int n )
{
	int i;
	for( i = 0; i < n; i++ ) {
		x[i] = 0.0;
	}
	return;
}
/*-----------------------------------------------------------*
 * procedure   mvr2r:                                        *
 *             ~~~~~~                                        *
 *  Transfer the contents of the vector x[] (real format)    *
 *  to the vector y[] (real format)                          *
 *-----------------------------------------------------------*/
void mvr2r(
    float x[], /* input : input vector  */
    float y[], /* output: output vector */
    int   n    /* input : vector size   */
)
{
	int i;
	for( i = 0; i < n; i++ ) {
		y[i] = x[i];
	}
	return;
}
/*-----------------------------------------------------------*
 * procedure   mvs2s:                                        *
 *             ~~~~~~                                        *
 *  Transfer the contents of the vector x[] (short format)   *
 *  to the vector y[] (short format)                         *
 *-----------------------------------------------------------*/
void mvs2s(
    short x[], /* input : input vector  */
    short y[], /* output: output vector */
    int   n    /* input : vector size   */
)
{
	int i;
	for( i = 0; i < n; i++ ) {
		y[i] = x[i];
	}
	return;
}
/*-----------------------------------------------------------*
 * procedure   mvs2s:                                        *
 *             ~~~~~~                                        *
 *  Transfer the contents of the vector x[] (short format)   *
 *  to the vector y[] (short format)                         *
 *-----------------------------------------------------------*/
void mvi2i(
    int x[], /* input : input vector  */
    int y[], /* output: output vector */
    int n    /* input : vector size   */
)
{
	int i;
	for( i = 0; i < n; i++ ) {
		y[i] = x[i];
	}
	return;
}
/*-----------------------------------------------------------*
 * procedure   mvr2s:                                        *
 *             ~~~~~~                                        *
 *  Transfer the contents of the vector x[] (real format)    *
 *  to the vector y[] (short format)                         *
 *-----------------------------------------------------------*/
void mvr2s(
    float x[], /* input : input vector  */
    short y[], /* output: output vector */
    int   n    /* input : vector size   */
)
{
	int   i;
	float temp;
	for( i = 0; i < n; i++ ) {
		temp = x[i];
		temp = (float)floor( temp + 0.5 );
		if( temp > 32767.0 )
			temp = 32767.0;
		if( temp < -32768.0 )
			temp = -32768.0;
		y[i] = (short)temp;
	}
	return;
}
/*-----------------------------------------------------------*
 * procedure   mvs2r:                                        *
 *             ~~~~~~                                        *
 *  Transfer the contents of the vector x[] (short format)   *
 *  to the vector y[] (reel format)                          *
 *-----------------------------------------------------------*/
void mvs2r(
    short x[], /* input : input vector  */
    float y[], /* output: output vector */
    int   n    /* input : vector size   */
)
{
	int i;
	for( i = 0; i < n; i++ ) {
		y[i] = (float)x[i];
	}
	return;
}

int get_nb_bits( short extension, short mode, short st_mode )
{
	int nb_bits;

	if( mode != 14 && mode != 15 ) /*prevent reading outside NBITS_CORE_AMR_WB buffer */
	{
		if( extension != 0 ) {
			nb_bits = NBITS_CORE[mode] + NBITS_BWE;
			if( st_mode >= 0 ) {
				nb_bits += ( StereoNbits[st_mode] + NBITS_BWE );
			}
		}
		else {
			nb_bits = NBITS_CORE_AMR_WB[mode];
		}
	}
	else {
		nb_bits = 0;
	}

	return nb_bits;
}
/*----------------------------------------------------------*
 * procedure pessimize:                                     *
 * A fake don-nothing routine to break program flow and     *
 * prevent optimisation                                     *
 *----------------------------------------------------------*/
void pessimize()
{
}
