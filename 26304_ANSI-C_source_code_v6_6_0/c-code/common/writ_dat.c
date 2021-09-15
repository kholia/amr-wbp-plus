/*----------------------------------------------------------------------
 *
 *  FUNCTION:   writ_data
 *
 *  PURPOSE:  round array of float data to 16-bit words and write to the
 *            file "fp"
 *
 *--------------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <stdio.h>
void writ_data(
    float data[], /* input : data              */
    int   size,   /* input : number of samples */
    FILE *fp      /* output: file pointer      */
)
{
	short data16[4 * L_FRAME_FSMAX];
	int   i;
	float temp;
	for( i = 0; i < size; i++ ) {
		temp = data[i];
		if( temp >= 0.0 )
			temp += 0.5;
		else
			temp -= 0.5;
		if( temp > 32767.0 )
			temp = 32767.0;
		if( temp < -32767.0 )
			temp = -32767.0;
		data16[i] = (short)temp;
	}
	fwrite( data16, sizeof( short ), size, fp );
	return;
}
