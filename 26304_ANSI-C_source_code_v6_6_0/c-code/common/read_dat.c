/*-----------------------------------------------------------------------
 * FUNCTION:   read_data
 *
 * PURPOSE:   reads  "size" 16-bit speech samples from the file "fp" and
 *            output to the array "data[]"
 *
 *-------------------------------------------------------------------------*/
#include "../include/amr_plus.h"
#include <stdio.h>
int read_data(    /* return: number of data successfully read */
    FILE *fp,     /* input : data file (16-bit words)         */
    float data[], /* output: speech data                      */
    int   size    /* input : number of samples                */
)
{
	int   i, n;
	short data16[4 * L_FRAME_FSMAX];
	n = fread( (void *)data16, sizeof( short ), size, fp );
	for( i = 0; i < n; i++ ) {
		data[i] = (float)data16[i];
	}
	for( i = n; i < size; i++ ) {
		data[i] = 0.0;
	}
	return n;
}
