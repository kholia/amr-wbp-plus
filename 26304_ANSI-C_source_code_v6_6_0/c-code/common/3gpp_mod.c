#include "../include/int3gpp.h"
#include "../lib_amr/enc_util.h"
#include "../lib_amr/typedef.h"
#include <math.h>
#include <memory.h>
#define L_FRAME16k 320    /* Frame size at 16kHz */
#define M16k 20           /* Order of LP filter */
#define L_WINDOW_PLUS 512 /* 448 low rate, 512 using EXTENSION_VA */
void E_LPC_isf_reorderPlus( float *isf, float min_dist, int n )
{
	int   i;
	float isf_min;
	isf_min = min_dist;
	for( i = 0; i < n - 1; i++ ) {
		if( isf[i] < isf_min ) {
			isf[i] = isf_min;
		}
		isf_min = isf[i] + min_dist;
	}
	return;
}
void E_UTIL_synthesisPlus( Float32 a[], Word32 m, Float32 x[], Float32 y[], Word32 l, Float32 mem[], Word32 update_m )
{
	Float32  buf[L_FRAME16k + M16k]; /* temporary synthesis buffer */
	Float32  s;
	Float32 *yy;
	Word32   i, j;
	/* copy initial filter states into synthesis buffer */
	memcpy( buf, mem, m * sizeof( Float32 ) );
	yy = &buf[m];
	for( i = 0; i < l; i++ ) {
		s = x[i];
		for( j = 1; j <= m; j++ ) {
			s -= a[j] * yy[i - j];
		}
		yy[i] = s;
		y[i] = s;
	}
	/* Update memory if required */
	if( update_m ) {
		memcpy( mem, &yy[l - m], m * sizeof( Float32 ) );
	}
	return;
}
void E_UTIL_residuPlus( Float32 *a, Word32 m, Float32 *x, Float32 *y, Word32 l )
{
	Float32 s;
	Word32  i, j;
	for( i = 0; i < l; i++ ) {
		s = x[i];
		for( j = 1; j <= m; j++ ) {
			s += a[j] * x[i - j];
		}
		y[i] = s;
	}
	return;
}
void E_UTIL_autocorrPlus( float *x, /* input : input signal */
    float *                      r, /* output: autocorrelations vector */
    int                          m, /* input : order of LP filter */
    int                          n, /* input : window size */
    float *                      fh /* input : analysis window */
)
{
	float  t[L_WINDOW_PLUS];
	float  s;
	Word16 i, j;
	for( i = 0; i < n; i++ ) {
		t[i] = x[i] * fh[i];
	}
	for( i = 0; i <= m; i++ ) {
		s = 0.0;
		for( j = 0; j < n - i; j++ ) {
			s += t[j] * t[j + i];
		}
		r[i] = s;
	}
	if( r[0] < 1.0 ) {
		r[0] = 1.0;
	}
	return;
}
