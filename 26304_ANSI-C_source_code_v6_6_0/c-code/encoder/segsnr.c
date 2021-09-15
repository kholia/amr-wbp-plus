/*_____________________________________________________________________
 |                                                                     |
 |  FUNCTION NAME segsnr                                               |
 |      Computes the segmential signal-to-noise ratio between the      |
 |      signal x and its estimate xe of length n samples. The segment  |
 |      length is nseg samples.                                        |
 |
 |  INPUT
 |      x[0:n-1]   Signal of n samples.
 |      xe[0:n-1]  Estimated signal of n samples.
 |      n          Signal length.
 |      nseg       Segment length, must be a submultiple of n.
 |
 |  RETURN VALUE
 |      snr        Segmential signal to noise ratio in dB.
 |_____________________________________________________________________|
*/
#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
float segsnr( float x[], float xe[], short n, short nseg )
{
	float snr = 0.0f;
	float signal, noise, error, fac;
	short i, j;
	for( i = 0; i < n; i += nseg ) {
		signal = 1e-6f;
		noise = 1e-6f;
		for( j = 0; j < nseg; j++ ) {
			signal += ( *x ) * ( *x );
			error = *x++ - *xe++;
			noise += error * error;
		}
		snr += (float)log10( (double)( signal / noise ) );
	}
	fac = ( (float)( 10 * nseg ) ) / (float)n;
	snr = fac * snr;
	if( snr < -99.0f ) {
		snr = -99.0f;
	}
	return ( snr );
}
