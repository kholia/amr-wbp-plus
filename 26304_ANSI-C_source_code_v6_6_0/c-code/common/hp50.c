#include "../include/amr_plus.h"
/*
 * hp50_12k8
 *
 * Function:
 *    2nd order high pass filter with nominal cut off frequency at 50 Hz.
 *
 * Returns:
 *    void
 */
void hp50_12k8( Float32 signal[], Word32 lg, Float32 mem[], Word32 fscale )
{
	Word16  i;
	Float32 x0, x1, x2, y0, y1, y2;
	Float32 a1, a2, b1, b2, frac;
	y1 = mem[0];
	y2 = mem[1];
	x0 = mem[2];
	x1 = mem[3];
	if( fscale >= FSCALE_DENOM ) {
		/* design in matlab (-3dB @ 20 Hz)
         [b,a] = butter(2, 20.0/6400.0, 'high');
         a = [1.00000000000000  -1.98611621154089   0.98621192916075];
         b = [0.99308203517541  -1.98616407035082   0.99308203517541];
         [b,a] = butter(2, 20.0/12800.0, 'high');
         a = [1.00000000000000  -1.99305802314321   0.99308203546221];
         b = [0.99653501465135  -1.99307002930271   0.99653501465135];
         difference:
         a = [               0  -0.00694181160232   0.00687010630146];
         b = [0.00345297947594  -0.00690595895189   0.00345297947594];
      */
		/* this approximation is enough (max overshoot = +0.05dB) */
		frac = ( ( Float32 )( fscale - FSCALE_DENOM ) ) / ( (Float32)FSCALE_DENOM );
		a1 = 1.98611621154089f + ( frac * 0.00694181160232f );
		a2 = -0.98621192916075f - ( frac * 0.00687010630146f );
		b1 = -1.98616407035082f - ( frac * 0.00690595895189f );
		b2 = 0.99308203517541f + ( frac * 0.00345297947594f );
		for( i = 0; i < lg; i++ ) {
			x2 = x1;
			x1 = x0;
			x0 = signal[i];
			y0 = ( y1 * a1 ) + ( y2 * a2 ) + ( x0 * b2 ) + ( x1 * b1 ) + ( x2 * b2 );
			signal[i] = y0;
			y2 = y1;
			y1 = y0;
		}
	}
	else /* -6dB 24Hz fs/2=6400Hz (3GPP AMR-WB coef) */
	{
		for( i = 0; i < lg; i++ ) {
			x2 = x1;
			x1 = x0;
			x0 = signal[i];
			y0 = y1 * 1.978881836F + y2 * -0.979125977F + x0 * 0.989501953F + x1 * -1.979003906F + x2 * 0.989501953F;
			signal[i] = y0;
			y2 = y1;
			y1 = y0;
		}
	}
	mem[0] = ( ( y1 > 1e-10 ) | ( y1 < -1e-10 ) ) ? y1 : 0;
	mem[1] = ( ( y2 > 1e-10 ) | ( y2 < -1e-10 ) ) ? y2 : 0;
	mem[2] = ( ( x0 > 1e-10 ) | ( x0 < -1e-10 ) ) ? x0 : 0;
	mem[3] = ( ( x1 > 1e-10 ) | ( x1 < -1e-10 ) ) ? x1 : 0;
	return;
}
