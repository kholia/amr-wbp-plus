#include "../include/amr_plus.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*_____________________________________________________________________
 |                                                                     |
 |                       Constant Definitions                          |
 |_____________________________________________________________________|
*/
#define ORDER_MAX 7
void fft3( float X[], float Y[], short n )
{
	float  Z[N_MAX / 6];
	float *Z0, *Z1, *Z2;
	float *z0, *z1, *z2;
	float *x;
	float *yre, *yim, *zre, *zim, *wre, *wim;
	short  m, step, sign, order;
	short  i, j, k;
	short  Ind[3];
	short *ind;
	/* Determine the order of the transform, the length of decimated  */
	/* transforms m, and the step for the sine and cosine tables.     */
	switch( n ) {
	case 48:
		order = 4;
		m = 16;
		step = 24;
		break;
	case 96:
		order = 5;
		m = 32;
		step = 12;
		break;
	case 192:
		order = 6;
		m = 64;
		step = 6;
		break;
	default:
		printf( " invalid fft3 size!\n" );
		exit( 0 );
	}
	/* Compose decimated sequences X[3i], X[3i+1],X[3i+2] */
	/* compute their FFT of length m.                                 */
	Z0 = &Z[0];
	z0 = &Z0[0];
	Z1 = &Z0[m];
	z1 = &Z1[0]; /* Z1 = &Z[ m];     */
	Z2 = &Z1[m];
	z2 = &Z2[0]; /* Z2 = &Z[2m];     */
	x = &X[0];
	for( i = 0; i < n / 3; i++ ) {
		*z0++ = *x++; /* Z0[i] = X[3i];   */
		*z1++ = *x++; /* Z1[i] = X[3i+1]; */
		*z2++ = *x++; /* Z2[i] = X[3i+2]; */
	}
	fft_rel( &Z0[0], m, order );
	fft_rel( &Z1[0], m, order );
	fft_rel( &Z2[0], m, order );
	/* Compute the DC coefficient and store it into Y[0]. Note that   */
	/* the variables Z0, ..., Z8, z0, ..., z8 are not needed after    */
	/* this.                                                          */
	*Y = *Z0 + *Z1 + *Z2;
	/* Initialize the index table, which points to the sine and       */
	/* cosine tables.                                                 */
	ind = &Ind[0];
	for( k = 1; k < 3; k++ ) {
		*ind++ = (short)( k * step );
	}
	/* Butterflies of order 3. */
	sign = 1;
	zre = &Z[1];
	zim = &Z[m - 1];
	yre = &Y[1];
	yim = &Y[n / 2 + 1];
	for( i = 0; i < 3; i++ ) {
		for( j = 1; j < m / 2; j++ ) {
			wre = &zre[0];
			wim = &zim[0];
			ind = &Ind[0];
			*yre = *wre;
			*yim = sign * ( *wim );
			for( k = 1; k < 3; k++ ) {
				wre += m;
				wim += m;
				*yre += ( *wre ) * t_cos[*ind] + sign * ( *wim ) * t_sin[*ind];
				*yim += -( *wre ) * t_sin[*ind] + sign * ( *wim ) * t_cos[*ind];
				*ind = (short)( *ind + k * step );
				if( *ind >= N_MAX ) {
					*ind -= N_MAX;
				}
				ind++;
			}
			yre++;
			zre += sign;
			yim++;
			zim -= sign;
		}
		wre = &zre[0];
		ind = &Ind[0];
		*yre = *wre;
		/*        *yim = 0.0; */
		if( i < 2 ) {
			*yim = 0.0;
		}
		for( k = 1; k < 3; k++ ) {
			wre += m;
			*yre += ( *wre ) * t_cos[*ind];
			/*            *yim += -(*wre)*t_sin[*ind]; */
			if( i < 2 ) {
				*yim += -( *wre ) * t_sin[*ind];
			}
			*ind = (short)( *ind + k * step );
			if( *ind >= N_MAX ) {
				*ind -= N_MAX;
			}
			ind++;
		}
		sign = (short)-sign;
		yre++;
		zre += sign;
		yim++;
		zim -= sign;
	}
	for( i = 0; i < n; i++ ) {
		Z[i] = Y[i];
	}
	Y[1] = Z[n / 2];
	for( i = 1; i < ( n / 2 ); i++ ) {
		Y[i * 2] = Z[i];
		Y[( i * 2 ) + 1] = Z[( n / 2 ) + i];
	}
	return;
}
void ifft3( float Y[], float X[], short n )
{
	short  i, j, k, m, step, order;
	short  Ind[3];
	short *ind;
	float  Z[N_MAX / 6];
	float *z, *zre, *zim;
	float *z0, *z1, *z2;
	float *yr0, *yr1;
	float *yi0, *yi1;
	float *yr0f, *yr1f;
	float *yi0f, *yi1f;
	float *yr0b;
	float *yi0b;
	float  scale;
	/* Determine the order of the transform, the length of decimated  */
	/* transforms m, and the step for the sine and cosine tables.     */
	switch( n ) {
	case 48:
		order = 4;
		m = 16;
		step = 24;
		break;
	case 96:
		order = 5;
		m = 32;
		step = 12;
		break;
	case 192:
		order = 6;
		m = 64;
		step = 6;
		break;
	default:
		printf( " invalid fft3 size!\n" );
		exit( 0 );
	}
	/* reordering : re[0], re[n/2], re[1], im[1], ..., re[n/2-1], im[n/2-1] */
	/* to           re[0], re[1], ... re[n/2],  im[1], im[2], ... im[n/2-1] */
	for( i = 0; i < n; i++ ) {
		Z[i] = Y[i];
	}
	Y[n / 2] = Z[1];
	for( i = 1; i < ( n / 2 ); i++ ) {
		Y[i] = Z[i * 2];
		Y[( n / 2 ) + i] = Z[( i * 2 ) + 1];
	}
	/* Initialize the fixed and the floating pointers. */
	yr0 = &Y[0];
	yr1 = &yr0[m]; /* = &Y[  m];     */
	yi0 = &Y[n / 2 + 1];
	yi1 = &Y[n / 2 + m];
	zre = &Z[0];
	zim = &Z[m - 1];
	yr0f = &yr0[0];
	yr0b = &yr1[-1];
	yr1f = &yr1[0];
	yi0f = &yi0[0];
	yi0b = &yi1[-1];
	yi1f = &yi0f[m];
	/* Compute the inverse butterflies. */
	/* p = 0*/
	*zre++ = *yr0f++ + 2 * ( *yr1f++ );
	for( i = 1; i < m / 2; i++ ) {
		*zre++ = *yr0f++ + *yr1f++ + *yr0b--;
		*zim-- = *yi0f++ + *yi1f++ - *yi0b--;
	}
	*zre = 2 * ( *yr0f ) + ( *yr1f );
	/* p=1,2 */
	for( k = 1; k < 3; k++ ) {
		ind = &Ind[0];
		*ind++ = 0;
		for( j = 1; j < 3; j++ ) {
			*ind = (short)( ind[-1] + m * step );
			if( *ind >= N_MAX ) {
				*ind -= N_MAX;
			}
			ind++;
		}
		z = &Z[k * m];
		*z = *yr0;
		ind = &Ind[1];
		*z += *yr1 * t_cos[*ind] - *yi1 * t_sin[*ind];
		ind++;
		*z += *yr1 * t_cos[*ind] + *yi1 * t_sin[*ind];
		zre = &z[1];
		zim = &z[m - 1];
		yr0f = &yr0[1];
		yi0f = &yi0[0];
		yr1f = &yr1[1];
		yi1f = &yi0f[m];
		yr0b = &yr1[-1];
		yi0b = &yi1[-1];
		for( i = 1; i < m / 2; i++ ) {
			ind = &Ind[0];
			for( j = 0; j < 3; j++ ) {
				*ind = (short)( *ind + step );
				if( *ind >= N_MAX ) {
					*ind -= N_MAX;
				}
				ind++;
			}
			ind = &Ind[0];
			*zre = *yr0f * t_cos[*ind] - *yi0f * t_sin[*ind];
			*zim = *yr0f * t_sin[*ind] + *yi0f * t_cos[*ind];
			ind++;
			*zre += *yr1f * t_cos[*ind] - *yi1f * t_sin[*ind];
			*zim += *yr1f * t_sin[*ind] + *yi1f * t_cos[*ind];
			ind++;
			*zre += *yr0b * t_cos[*ind] + *yi0b * t_sin[*ind];
			*zim += *yr0b * t_sin[*ind] - *yi0b * t_cos[*ind];
			yr0b--;
			yi0b--;
			yr0f++;
			yi0f++;
			yr1f++;
			yi1f++;
			zre++;
			zim--;
		}
		ind = &Ind[0];
		for( j = 0; j < 3; j++ ) {
			*ind = (short)( *ind + step );
			if( *ind >= N_MAX ) {
				*ind -= N_MAX;
			}
			ind++;
		}
		ind = &Ind[0];
		*zre = *yr0f * t_cos[*ind] - *yi0f * t_sin[*ind];
		ind++;
		*zre += *yr1f * ( t_cos[*ind] - t_sin[*ind] );
		ind++;
		*zre += *yr0b * t_cos[*ind] + *yi0b * t_sin[*ind];           /* I have some doubts here, must check !!!!*/
		                                                             /* Update the table step. */
		step = (short)( step + 3 * ( 1 << ( ORDER_MAX - order ) ) ); /* triple the step if we want to use va's cosine table*/
	}
	/* Compute the inverse FFT for all nine blocks. */
	z0 = &Z[0];
	z1 = &z0[m]; /* z1 = &Z[ m];     */
	z2 = &z1[m]; /* z2 = &Z[2m];     */
	ifft_rel( &z0[0], m, order );
	ifft_rel( &z1[0], m, order );
	ifft_rel( &z2[0], m, order );
	/* Decimation and scaling, scale = 1/9m. */
	scale = ( (float)1 / 3 ) / ( (float)m );
	for( i = 0; i < n / 3; i++ ) {
		*X++ = ( *z0++ ) * scale;
		*X++ = ( *z1++ ) * scale;
		*X++ = ( *z2++ ) * scale;
	}
}
