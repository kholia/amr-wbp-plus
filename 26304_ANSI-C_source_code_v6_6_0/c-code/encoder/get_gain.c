#include "../include/amr_plus.h"
float get_gain( /* output: codebook gain (adaptive or fixed) */
    float x[],  /* input : target signal                     */
    float y[],  /* input : filtered codebook excitation      */
    int   n     /* input : segment length                    */
)
{
	float corr = 0.0f, ener = 1e-6f;
	short i;
	for( i = 0; i < n; i++ ) {
		corr += x[i] * y[i];
		ener += y[i] * y[i];
	}
	return ( corr / ener );
}
