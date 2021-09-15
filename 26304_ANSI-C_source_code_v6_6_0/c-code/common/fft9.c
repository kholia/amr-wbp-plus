/*_____________________________________________________________________
 |                                                                     |
 |                        File Inclusions                              |
 |_____________________________________________________________________|
*/
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
/*_____________________________________________________________________
 |                                                                     |
 |  FUNCTION NAME fft9                                                 |
 |      Radix-9 FFT for real-valued sequences of length 576 or 1152.
 |      The  input  sequence  is  first decimated by nine,  and  the
 |      split-radix  algorithm described in [1]  are applied to  the
 |      decimated sequences.  The resulting nine separate transforms
 |      are then combined.
 |
 |      The function requires sine and cosine tables t_sin and t_cos,
 |      and constants N_MAX = 1152 and ORDER_MAX = log2(N_MAX/9). The
 |      table entries are  defined as sin(2*pi*i) and cos(2*pi*i) for
 |      i = 0, 1, ..., N_MAX.
 |
 |  INPUT
 |      X[0:n-1]  Input sequence.
 |      n         Number of samples in the sequence, must be 288,
 |                576, or 1152.
 |
 |  OUTPUT
 |      Y[0:n-1]  Transform coeffients in the order re[0], re[n/2],
 |                re[1], im[1], re[2], im[2], ... re[n/2-1], im[n/2-1].
 |_____________________________________________________________________|
*/
void fft9( float X[], float Y[], short n )
{
	float  Z[N_MAX];
	float *Z0, *Z1, *Z2, *Z3, *Z4, *Z5, *Z6, *Z7, *Z8;
	float *z0, *z1, *z2, *z3, *z4, *z5, *z6, *z7, *z8;
	float *x;
	float *yre, *yim, *zre, *zim, *wre, *wim;
	short  m, step, sign, order;
	short  i, j, k;
	short  Ind[9];
	short *ind;
	short  temp;
	/* Determine the order of the transform, the length of decimated  */
	/* transforms m, and the step for the sine and cosine tables.     */
	switch( n ) {
	case 72:
		order = 3;
		m = 8;
		step = 16;
		break;
	case 144:
		order = 4;
		m = 16;
		step = 8;
		break;
	case 288:
		order = 5;
		m = 32;
		step = 4;
		break;
	case 576:
		order = 6;
		m = 64;
		step = 2;
		break;
	case 1152:
		order = 7;
		m = 128;
		step = 1;
		break;
	default:
		printf( " invalid fft9 size!\n" );
		exit( 0 );
	}
	/* Compose decimated sequences X[9i], X[9i+1], ..., X[9i+4] and   */
	/* compute their FFT of length m.                                 */
	Z0 = &Z[0];
	z0 = &Z0[0];
	Z1 = &Z0[m];
	z1 = &Z1[0]; /* Z1 = &Z[ m];     */
	Z2 = &Z1[m];
	z2 = &Z2[0]; /* Z2 = &Z[2m];     */
	Z3 = &Z2[m];
	z3 = &Z3[0]; /* Z3 = &Z[3m];     */
	Z4 = &Z3[m];
	z4 = &Z4[0]; /* Z4 = &Z[4m];     */
	Z5 = &Z4[m];
	z5 = &Z5[0]; /* Z5 = &Z[5m];     */
	Z6 = &Z5[m];
	z6 = &Z6[0]; /* Z6 = &Z[6m];     */
	Z7 = &Z6[m];
	z7 = &Z7[0]; /* Z7 = &Z[7m];     */
	Z8 = &Z7[m];
	z8 = &Z8[0]; /* Z8 = &Z[8m];     */
	x = &X[0];
	for( i = 0; i < n / 9; i++ ) {
		*z0++ = *x++; /* Z0[i] = X[9i];   */
		*z1++ = *x++; /* Z1[i] = X[9i+1]; */
		*z2++ = *x++; /* Z2[i] = X[9i+2]; */
		*z3++ = *x++; /* Z3[i] = X[9i+3]; */
		*z4++ = *x++; /* Z4[i] = X[9i+4]; */
		*z5++ = *x++; /* Z5[i] = X[9i+5]; */
		*z6++ = *x++; /* Z6[i] = X[9i+6]; */
		*z7++ = *x++; /* Z7[i] = X[9i+7]; */
		*z8++ = *x++; /* Z8[i] = X[9i+8]; */
	}
	fft_rel( &Z0[0], m, order );
	fft_rel( &Z1[0], m, order );
	fft_rel( &Z2[0], m, order );
	fft_rel( &Z3[0], m, order );
	fft_rel( &Z4[0], m, order );
	fft_rel( &Z5[0], m, order );
	fft_rel( &Z6[0], m, order );
	fft_rel( &Z7[0], m, order );
	fft_rel( &Z8[0], m, order );
	/* Compute the DC coefficient and store it into Y[0]. Note that   */
	/* the variables Z0, ..., Z8, z0, ..., z8 are not needed after    */
	/* this.                                                          */
	*Y = *Z0 + *Z1 + *Z2 + *Z3 + *Z4 + *Z5 + *Z6 + *Z7 + *Z8;
	/* Initialize the index table, which points to the sine and       */
	/* cosine tables.                                                 */
	ind = &Ind[0];
	for( k = 1; k < 9; k++ ) {
		*ind++ = (short)( k * step );
	}
	/* Butterflies of order 9. */
	/* EXAMPLE RADIX5:                                                */
	/* ~~~~~~~~~~~~~~~                                                */
	/* Transform coefficients in Y are computed so that the pointer   */
	/* yre goes over Y[1:n/2] = Y[1:5m],  and the pointer yim  over   */
	/* Y[n/2+1:n-1] = Y[5m+1:10m-1].  The pointers  zre and zim run   */
	/* over the butterflies in Z according to the following table.    */
	/*                                                                */
	/* ===================================================            */
	/*   i     yre           yim       zre         zim                */
	/* ===================================================            */
	/*   0   Y[   1: m]  Y[10m-1:9m]  Z[1:m]    Z[2m-1:m]             */
	/*   1   Y[ m+1:2m]  Y[ 9m-1:8m]  Z[m-1:0]  Z[m+1:2m]             */
	/*   2   Y[2m+1:3m]  Y[ 8m-1:7m]  Z[1:m]    Z[2m-1:m]             */
	/*   3   Y[3m+1:4m]  Y[ 7m-1:6m]  Z[m-1:0]  Z[m+1:2m]             */
	/*   4   Y[4m+1:5m]  Y[ 6m-1:5m]  Z[1:m]    Z[2m-1:m]             */
	/* ===================================================            */
	sign = 1;
	zre = &Z[1];
	zim = &Z[m - 1];
	yre = &Y[1];
	yim = &Y[n / 2 + 1];
	for( i = 0; i < 9; i++ ) {
		for( j = 1; j < m / 2; j++ ) {
			wre = &zre[0];
			wim = &zim[0];
			ind = &Ind[0];
			*yre = *wre;
			*yim = sign * ( *wim );
			for( k = 1; k < 9; k++ ) {
				wre += m;
				wim += m;
				temp = *ind;
				*yre += ( *wre ) * t_cos[temp] + sign * ( *wim ) * t_sin[temp];
				*yim += -( *wre ) * t_sin[temp] + sign * ( *wim ) * t_cos[temp];
				temp = (short)( temp + k * step );
				if( temp >= N_MAX ) {
					temp -= N_MAX;
				}
				*ind = temp;
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
		if( i < 8 ) {
			*yim = 0.0;
		}
		for( k = 1; k < 9; k++ ) {
			wre += m;
			*yre += ( *wre ) * t_cos[*ind];
			/*            *yim += -(*wre)*t_sin[*ind]; */
			if( i < 8 ) {
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
	/* reordering : re[0], re[1], ... re[n/2],  im[1], im[2], ... im[n/2-1] */
	/* to           re[0], re[n/2], re[1], im[1], ..., re[n/2-1], im[n/2-1] */
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
/*_____________________________________________________________________
 |                                                                     |
 |  FUNCTION NAME fft_rel                                              |
 |      Computes the split-radix FFT in place for the real-valued
 |      signal x of length n.  The algorithm has been ported from
 |      the Fortran code of [1].
 |
 |      The function  needs sine and cosine tables t_sin and t_cos,
 |      and the constant N_MAX.  The table  entries  are defined as
 |      sin(2*pi*i) and cos(2*pi*i) for i = 0, 1, ..., N_MAX-1. The
 |      implementation  assumes  that any entry  will not be needed
 |      outside the tables. Therefore, N_MAX and n must be properly
 |      set.  The function has been tested  with the values n = 288,
 |      576 and N_MAX = 1152.
 |
 |      References
 |      [1] H.V. Sorensen,  D.L. Jones, M.T. Heideman, C.S. Burrus,
 |          "Real-valued fast  Fourier transform  algorithm,"  IEEE
 |          Trans. on Signal Processing,  Vol.35, No.6, pp 849-863,
 |          1987.
 |
 |  INPUT
 |      x[0:n-1]  Input sequence.
 |      n         Number of samples in the sequence.
 |      m         m = log2(n).
 |
 |  OUTPUT
 |      x[0:n-1]  Transform coeffients in the order re[0], re[1],
 |                ..., re[n/2], im[n/2-1], ..., im[1].
 |_____________________________________________________________________|
*/
void fft_rel( float x[], short n, short m )
{
	short  i, j, k, n1, n2, n4;
	short  i1, i2, i3, i4;
	short  step, ind;
	float  xt, t1, t2;
	float *x0, *x1;
	/* Digit reverse counter. */
	j = 0;
	x0 = &x[0];
	for( i = 0; i < n - 1; i++ ) {
		if( i < j ) {
			xt = x[j];  /*   xt = x[j] */
			x[j] = *x0; /* x[j] = x[i] */
			*x0 = xt;   /* x[i] = xt   */
		}
		x0++;
		k = (short)( n / 2 );
		while( k <= j ) {
			j = (short)( j - k );
			k = (short)( k >> 1 );
			/* k<=j comparison */
		}
		j = (short)( j + k );
	}
	/* Length two butterflies. */
	x0 = &x[0];
	x1 = &x[1];
	for( i = 0; i < n / 2; i++ ) {
		xt = *x0;       /* xt      = x[2i]        */
		*x0 = xt + *x1; /* x[2i]   = xt - x[2i+1] */
		*x1 = xt - *x1; /* x[2i+1] = xt - x[2i+1] */
		x0++;
		x0++;
		x1++;
		x1++;
	}
	/* Other butterflies. */
	/* The implementation described in [1] has been changed by using  */
	/* table lookup for evaluating sine and cosine functions.  The    */
	/* variable ind and its increment step are needed to access table */
	/* entries.  Note that this implementation assumes n4 to be so    */
	/* small that ind will never exceed the table.  Thus the input    */
	/* argument n and the constant N_MAX must be set properly.        */
	n2 = 1;
	for( k = 2; k <= m; k++ ) {
		n4 = n2;
		n2 = (short)( n4 << 1 );
		n1 = (short)( n2 << 1 );
		step = (short)( N_MAX / n1 );
		for( i = 0; i < n; i = (short)( i + n1 ) ) {
			xt = x[i];
			x[i] = xt + x[i + n2];
			x[i + n2] = xt - x[i + n2];
			x[i + n2 + n4] = -x[i + n2 + n4];
			ind = step;
			for( j = 1; j < n4; j++ ) {
				short temp;
				i1 = (short)( i + j );
				temp = (short)( i - j );
				i2 = (short)( temp + n2 );
				i3 = (short)( i1 + n2 );
				i4 = (short)( temp + n1 );
				t1 = x[i3] * t_cos[ind] + x[i4] * t_sin[ind];
				t2 = x[i3] * t_sin[ind] - x[i4] * t_cos[ind];
				x[i4] = x[i2] - t2;
				x[i3] = -x[i2] - t2;
				x[i2] = x[i1] - t1;
				x[i1] = x[i1] + t1;
				ind = (short)( ind + step );
			}
		}
	}
}
/*_____________________________________________________________________
 |                                                                     |
 |  FUNCTION NAME ifft9                                                |
 |      Inverse Radix-9 FFT for real-valued sequences of length 288,   |
 |      576 or 1152. The function computes first inverse butterflies   |
 |      for obtaining transform coefficients of sequencies decimated
 |      by 9.  The inverse split-radix FFT [1] is applied separately
 |      to these nine coefficient blocks and the outcomes are
 |      combined.
 |
 |      The function requires sine and cosine tables t_sin and t_cos,
 |      and constants N_MAX = 1152 and ORDER_MAX = log2(N_MAX/9). The
 |      table entries are  defined as sin(2*pi*i) and cos(2*pi*i) for
 |      i = 0, 1, ..., N_MAX-1.
 |
 |  INPUT
 |      Y[0:n-1]  Transform coeffients in the order re[0], re[n/2],
 |                re[1], im[1], re[2], im[2], ... re[n/2-1], im[n/2-1].
 |      n         Number of transform coefficients, must be 288, 576,
 |                or 1152.
 |
 |  OUTPUT
 |      X[0:n-1]  Output sequence.
 |_____________________________________________________________________|
*/
void ifft9( float Y[], float X[], short n )
{
	short  i, j, k, m, step, order;
	short  Ind[9];
	short *ind;
	float  Z[N_MAX];
	float *z, *zre, *zim;
	float *z0, *z1, *z2, *z3, *z4, *z5, *z6, *z7, *z8;
	float *yr0, *yr1, *yr2, *yr3, *yr4;
	float *yi0, *yi1, *yi2, *yi3, *yi4;
	float *yr0f, *yr1f, *yr2f, *yr3f, *yr4f;
	float *yi0f, *yi1f, *yi2f, *yi3f, *yi4f;
	float *yr0b, *yr1b, *yr2b, *yr3b;
	float *yi0b, *yi1b, *yi2b, *yi3b;
	float  scale;
	/* Determine the order of the transform, the length of decimated  */
	/* transforms m, and the step for the sine and cosine tables.     */
	switch( n ) {
	case 72:
		order = 3;
		m = 8;
		step = 16;
		break;
	case 144:
		order = 4;
		m = 16;
		step = 8;
		break;
	case 288:
		order = 5;
		m = 32;
		step = 4;
		break;
	case 576:
		order = 6;
		m = 64;
		step = 2;
		break;
	case 1152:
		order = 7;
		m = 128;
		step = 1;
		break;
	default:
		printf( " invalid ifft9 size!\n" );
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
	/* EXAMPLE RADIX5:                                                */
	/* ~~~~~~~~~~~~~~~                                                */
	/* The following table depicts indexing and locations of pointers */
	/* in an illustrative example case as n = 20 and m = n/5 = 4.     */
	/* The pointers yr0, yr1, yr2, yi0, yi1, yi2 are anchored to the  */
	/* beginning of the blocks in the coefficient vector. They will   */
	/* not be changed during the computation. The coefficient vector  */
	/* and the fixed pointers are shown in the second and third       */
	/* column of the table. The floating pointers yr0f, yr1f, yr2f,   */
	/* yi0f, yi1f, yi2f go forward over the corresponding blocks.     */
	/* Correspondingly yr0b, yr1b, yi0b, yi1b run backwards. These    */
	/* pointers will be repositioned during the algorithm.            */
	/*                                                                */
	/*   ==================================                           */
	/*   index  coeff    block     pointers                           */
	/*   ==================================                           */
	/*     0    re[0]    &yr0[0]    yr0f->                            */
	/*     1    re[1]                                                 */
	/*     2    re[2]                                                 */
	/*     3    re[3]               <-yr0b                            */
	/*     4    re[4]    &yr1[0]    yr1f->                            */
	/*     5    re[5]                                                 */
	/*     6    re[6]                                                 */
	/*     7    re[7]               <-yr1b                            */
	/*     8    re[8]    &yr2[0]    yr2f->                            */
	/*     9    re[9]                                                 */
	/*     10   re[10]                                                */
	/*     11   im[1]    &yi0[0]    yi0f->                            */
	/*     12   im[2]                                                 */
	/*     13   im[3]               <-yi0b                            */
	/*     14   im[4]    &yi1[0]                                      */
	/*     15   im[5]               yi1f->                            */
	/*     16   im[6]                                                 */
	/*     17   im[7]               <-yi1b                            */
	/*     18   im[8]    &yi2[0]                                      */
	/*     19   im[9]               yi2f->                            */
	/*   ==================================                           */
	/* Initialize the fixed and the floating pointers. */
	yr0 = &Y[0];
	yr1 = &yr0[m]; /* = &Y[  m];     */
	yr2 = &yr1[m]; /* = &Y[2*m];     */
	yr3 = &yr2[m];
	yr4 = &yr3[m];
	yi0 = &Y[n / 2 + 1];
	yi1 = &Y[n / 2 + m];
	yi2 = &Y[n / 2 + 2 * m];
	yi3 = &Y[n / 2 + 3 * m];
	yi4 = &Y[n / 2 + 4 * m];
	zre = &Z[0];
	zim = &Z[m - 1];
	yr0f = &yr0[0];
	yr0b = &yr1[-1];
	yr1f = &yr1[0];
	yr1b = &yr2[-1];
	yr2f = &yr2[0];
	yr2b = &yr3[-1];
	yr3f = &yr3[0];
	yr3b = &yr4[-1];
	yr4f = &yr4[0];
	yi0f = &yi0[0];
	yi0b = &yi1[-1];
	yi1f = &yi0f[m];
	yi1b = &yi2[-1];
	yi2f = &yi1f[m];
	yi2b = &yi3[-1];
	yi3f = &yi2f[m];
	yi3b = &yi4[-1];
	yi4f = &yi3f[m];
	/* Compute the inverse butterflies. */
	*zre++ = *yr0f++ + 2 * ( *yr1f++ ) + 2 * ( *yr2f++ ) + 2 * ( *yr3f++ ) + 2 * ( *yr4f++ );
	for( i = 1; i < m / 2; i++ ) {
		*zre++ = *yr0f++ + *yr1f++ + *yr2f++ + *yr3f++ + *yr4f++ + *yr0b-- + *yr1b-- + *yr2b-- + *yr3b--;
		*zim-- = *yi0f++ + *yi1f++ + *yi2f++ + *yi3f++ + *yi4f++ - *yi0b-- - *yi1b-- - *yi2b-- - *yi3b--;
	}
	*zre = 2 * ( *yr0f ) + 2 * ( *yr1f ) + 2 * ( *yr2f ) + 2 * ( *yr3f ) + *yr4f;
	for( k = 1; k < 9; k++ ) {
		ind = &Ind[0];
		*ind++ = 0;
		for( j = 1; j < 9; j++ ) {
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
		*z += *yr2 * t_cos[*ind] - *yi2 * t_sin[*ind];
		ind++;
		*z += *yr3 * t_cos[*ind] - *yi3 * t_sin[*ind];
		ind++;
		*z += *yr4 * t_cos[*ind] - *yi4 * t_sin[*ind];
		ind++;
		*z += *yr4 * t_cos[*ind] + *yi4 * t_sin[*ind];
		ind++;
		*z += *yr3 * t_cos[*ind] + *yi3 * t_sin[*ind];
		ind++;
		*z += *yr2 * t_cos[*ind] + *yi2 * t_sin[*ind];
		ind++;
		*z += *yr1 * t_cos[*ind] + *yi1 * t_sin[*ind];
		zre = &z[1];
		zim = &z[m - 1];
		yr0f = &yr0[1];
		yi0f = &yi0[0];
		yr1f = &yr1[1];
		yi1f = &yi0f[m];
		yr2f = &yr2[1];
		yi2f = &yi1f[m];
		yr3f = &yr3[1];
		yi3f = &yi2f[m];
		yr4f = &yr4[1];
		yi4f = &yi3f[m];
		yr3b = &yr4[-1];
		yi3b = &yi4[-1];
		yr2b = &yr3[-1];
		yi2b = &yi3[-1];
		yr1b = &yr2[-1];
		yi1b = &yi2[-1];
		yr0b = &yr1[-1];
		yi0b = &yi1[-1];
		for( i = 1; i < m / 2; i++ ) {
			ind = &Ind[0];
			for( j = 0; j < 9; j++ ) {
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
			*zre += *yr2f * t_cos[*ind] - *yi2f * t_sin[*ind];
			*zim += *yr2f * t_sin[*ind] + *yi2f * t_cos[*ind];
			ind++;
			*zre += *yr3f * t_cos[*ind] - *yi3f * t_sin[*ind];
			*zim += *yr3f * t_sin[*ind] + *yi3f * t_cos[*ind];
			ind++;
			*zre += *yr4f * t_cos[*ind] - *yi4f * t_sin[*ind];
			*zim += *yr4f * t_sin[*ind] + *yi4f * t_cos[*ind];
			ind++;
			yr0f++;
			yi0f++;
			yr1f++;
			yi1f++;
			yr2f++;
			yi2f++;
			yr3f++;
			yi3f++;
			yr4f++;
			yi4f++;
			*zre += *yr3b * t_cos[*ind] + *yi3b * t_sin[*ind];
			*zim += *yr3b * t_sin[*ind] - *yi3b * t_cos[*ind];
			ind++;
			*zre += *yr2b * t_cos[*ind] + *yi2b * t_sin[*ind];
			*zim += *yr2b * t_sin[*ind] - *yi2b * t_cos[*ind];
			ind++;
			*zre += *yr1b * t_cos[*ind] + *yi1b * t_sin[*ind];
			*zim += *yr1b * t_sin[*ind] - *yi1b * t_cos[*ind];
			ind++;
			*zre += *yr0b * t_cos[*ind] + *yi0b * t_sin[*ind];
			*zim += *yr0b * t_sin[*ind] - *yi0b * t_cos[*ind];
			yr0b--;
			yi0b--;
			yr1b--;
			yi1b--;
			yr2b--;
			yi2b--;
			yr3b--;
			yi3b--;
			zre++;
			zim--;
		}
		ind = &Ind[0];
		for( j = 0; j < 9; j++ ) {
			*ind = (short)( *ind + step );
			if( *ind >= N_MAX ) {
				*ind -= N_MAX;
			}
			ind++;
		}
		ind = &Ind[0];
		*zre = *yr0f * t_cos[*ind] - *yi0f * t_sin[*ind];
		ind++;
		*zre += *yr1f * t_cos[*ind] - *yi1f * t_sin[*ind];
		ind++;
		*zre += *yr2f * t_cos[*ind] - *yi2f * t_sin[*ind];
		ind++;
		*zre += *yr3f * t_cos[*ind] - *yi3f * t_sin[*ind];
		ind++;
		*zre += *yr4f * ( t_cos[*ind] - t_sin[*ind] );
		ind++;
		*zre += *yr3b * t_cos[*ind] + *yi3b * t_sin[*ind];
		ind++;
		*zre += *yr2b * t_cos[*ind] + *yi2b * t_sin[*ind];
		ind++;
		*zre += *yr1b * t_cos[*ind] + *yi1b * t_sin[*ind];
		ind++;
		*zre += *yr0b * t_cos[*ind] + *yi0b * t_sin[*ind];
		/* Update the table step. */
		step = (short)( step + ( 1 << ( ORDER_MAX - order ) ) );
	}
	/* Compute the inverse FFT for all nine blocks. */
	z0 = &Z[0];
	z1 = &z0[m]; /* z1 = &Z[ m];     */
	z2 = &z1[m]; /* z2 = &Z[2m];     */
	z3 = &z2[m]; /* z3 = &Z[3m];     */
	z4 = &z3[m]; /* z4 = &Z[4m];     */
	z5 = &z4[m]; /* z5 = &Z[5m];     */
	z6 = &z5[m]; /* z6 = &Z[6m];     */
	z7 = &z6[m]; /* z7 = &Z[7m];     */
	z8 = &z7[m]; /* z8 = &Z[8m];     */
	ifft_rel( &z0[0], m, order );
	ifft_rel( &z1[0], m, order );
	ifft_rel( &z2[0], m, order );
	ifft_rel( &z3[0], m, order );
	ifft_rel( &z4[0], m, order );
	ifft_rel( &z5[0], m, order );
	ifft_rel( &z6[0], m, order );
	ifft_rel( &z7[0], m, order );
	ifft_rel( &z8[0], m, order );
	/* Decimation and scaling, scale = 1/9m. */
	scale = ( (float)1 / 9 ) / ( (float)m );
	for( i = 0; i < n / 9; i++ ) {
		*X++ = ( *z0++ ) * scale;
		*X++ = ( *z1++ ) * scale;
		*X++ = ( *z2++ ) * scale;
		*X++ = ( *z3++ ) * scale;
		*X++ = ( *z4++ ) * scale;
		*X++ = ( *z5++ ) * scale;
		*X++ = ( *z6++ ) * scale;
		*X++ = ( *z7++ ) * scale;
		*X++ = ( *z8++ ) * scale;
	}
}
/*_____________________________________________________________________
 |                                                                     |
 |  FUNCTION NAME ifft_rel                                             |
 |      Computes the inverse split-radix FFT in place for the real-
 |      valued signal x of length n.  The algorithm has been ported
 |      from the Fortran code presented in [1].
 |
 |      The function  needs sine and cosine tables t_sin and t_cos,
 |      and the constant N_MAX.  The table  entries  are defined as
 |      sin(2*pi*i) and cos(2*pi*i) for i = 0, 1, ..., N_MAX-1. The
 |      implementation  assumes  that any entry  will not be needed
 |      outside the tables. Therefore, N_MAX and n must be properly
 |      set.  The function has been  tested with the values n = 16,
 |      32, 64, 128, 256, and N_MAX = 1280.
 |
 |      References
 |      [1] H.V. Sorensen,  D.L. Jones, M.T. Heideman, C.S. Burrus,
 |          "Real-valued fast  Fourier transform  algorithm,"  IEEE
 |          Trans. on Signal Processing,  Vol.35, No.6, pp 849-863,
 |          1987.
 |
 |  INPUT
 |      x[0:n-1]  Transform coeffients in the order re[0], re[1],
 |                ..., re[n/2], im[n/2-1], ..., im[1].
 |      n         Number of transform coefficients.
 |      m         m = log2(n).
 |
 |  OUTPUT
 |      x[0:n-1]  Output sequence.
 |_____________________________________________________________________|
*/
#define SQRT2 1.414213562373095048801688724209698078569671f
void ifft_rel( float x[], short n, short m )
{
	short  i, j, k;
	short  i0, i1, i2, i3, i4, i5, i6, i7, i8;
	short  n1, n2, n4, n8;
	short  is, id;
	short  step, ind;
	float *x0;
	float  xt;
	float  r1;
	float  t1, t2, t3, t4, t5;
	float  cc1, cc3, ss1, ss3;
	/* This patch is needed for converting the Fortran indexing to    */
	/* the C indexing.                                                */
	x = &x[-1];
	/* L-shaped butterflies.   */
	n2 = (short)( 2 * n );
	for( k = 1; k < m; k++ ) {
		is = 0;
		id = n2;
		n2 = (short)( n2 >> 1 );
		n4 = (short)( n2 >> 2 );
		n8 = (short)( n4 >> 1 );
		while( is < n - 1 ) {
			for( i = is; i < n; i = (short)( i + id ) ) {
				i1 = (short)( i + 1 );
				i2 = (short)( i1 + n4 );
				i3 = (short)( i2 + n4 );
				i4 = (short)( i3 + n4 );
				t1 = x[i1] - x[i3];
				x[i1] = x[i1] + x[i3];
				x[i2] = 2.0f * x[i2];
				x[i3] = t1 - 2.0f * x[i4];
				x[i4] = t1 + 2.0f * x[i4];
				if( n4 != 1 ) {
					i1 = (short)( i1 + n8 );
					i2 = (short)( i2 + n8 );
					i3 = (short)( i3 + n8 );
					i4 = (short)( i4 + n8 );
					t1 = x[i2] - x[i1];
					t2 = x[i4] + x[i3];
					x[i1] = x[i1] + x[i2];
					x[i2] = x[i4] - x[i3];
					x[i3] = SQRT2 * ( -t2 - t1 );
					x[i4] = SQRT2 * ( -t2 + t1 );
				}
			}
			is = (short)( 2 * id - n2 );
			id = (short)( 4 * id );
		}
		step = (short)( N_MAX / n2 );
		ind = step;
		for( j = 2; j <= n8; j++ ) {
			cc1 = t_cos[ind];
			ss1 = t_sin[ind];
			cc3 = t_cos[3 * ind];
			ss3 = t_sin[3 * ind];
			ind = (short)( ind + step );
			is = 0;
			id = (short)( 2 * n2 );
			while( is < n - 1 ) {
				for( i = is; i < n; i = (short)( i + id ) ) {
					i1 = (short)( i + j );
					i2 = (short)( i1 + n4 );
					i3 = (short)( i2 + n4 );
					i4 = (short)( i3 + n4 );
					i5 = (short)( i + n4 - j + 2 );
					i6 = (short)( i5 + n4 );
					i7 = (short)( i6 + n4 );
					i8 = (short)( i7 + n4 );
					t1 = x[i1] - x[i6];
					x[i1] = x[i1] + x[i6];
					t2 = x[i5] - x[i2];
					x[i5] = x[i2] + x[i5];
					t3 = x[i8] + x[i3];
					x[i6] = x[i8] - x[i3];
					t4 = x[i4] + x[i7];
					x[i2] = x[i4] - x[i7];
					t5 = t1 - t4;
					t1 = t1 + t4;
					t4 = t2 - t3;
					t2 = t2 + t3;
					x[i3] = t5 * cc1 + t4 * ss1;
					x[i7] = -t4 * cc1 + t5 * ss1;
					x[i4] = t1 * cc3 - t2 * ss3;
					x[i8] = t2 * cc3 + t1 * ss3;
				}
				is = (short)( 2 * id - n2 );
				id = (short)( 4 * id );
			}
		}
	}
	/* Length two butterflies. */
	is = 1;
	id = 4;
	while( is < n ) {
		for( i0 = is; i0 <= n; i0 = (short)( i0 + id ) ) {
			i1 = (short)( i0 + 1 );
			r1 = x[i0];
			x[i0] = r1 + x[i1];
			x[i1] = r1 - x[i1];
			/* two ptrs: x[i0] and x[i1], incrment by id, not 1 */
		}
		is = (short)( 2 * id - 1 );
		id = (short)( 4 * id );
	}
	/* Digit reverse counter. */
	j = 1;
	n1 = (short)( n - 1 );
	x0 = &x[1];
	for( i = 1; i < n; i++ ) {
		if( i < j ) {
			xt = x[j];  /* xt   = x[j] */
			x[j] = *x0; /* x[j] = x[i] */
			*x0 = xt;   /* x[i] = xt   */
		}
		x0++;
		k = (short)( n >> 1 );
		while( k < j ) {
			j = (short)( j - k );
			k = (short)( k >> 1 );
		}
		j = (short)( j + k );
	}
}
