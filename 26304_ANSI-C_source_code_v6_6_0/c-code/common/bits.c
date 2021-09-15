#include "../include/amr_plus.h"
#define MASK 0x0001
/*---------------------------------------------------------------------------*
 * function:  bin2int                                                        *
 *            ~~~~~~~                                                        *
 * Read "no_of_bits" bits from the array bitstream[] and convert to integer  *
 *--------------------------------------------------------------------------*/
int bin2int(           /* output: recovered integer value              */
    int    no_of_bits, /* input : number of bits associated with value */
    short *bitstream   /* input : address where bits are read          */
)
{
	int value, i;
	value = 0;
	for( i = 0; i < no_of_bits; i++ ) {
		value <<= 1;
		value += (int)( ( *bitstream++ ) & MASK );
	}
	return ( value );
}
/*---------------------------------------------------------------------------*
 * function:  int2bin                                                        *
 *            ~~~~~~~                                                        *
 * Convert integer to binary and write the bits to the array bitstream[].    *
 * Most significant bits (MSB) are output first                              *
 *--------------------------------------------------------------------------*/
void int2bin(
    int    value,      /* input : value to be converted to binary      */
    int    no_of_bits, /* input : number of bits associated with value */
    short *bitstream   /* output: address where bits are written       */
)
{
	short *pt_bitstream;
	int    i;
	pt_bitstream = bitstream + no_of_bits;
	for( i = 0; i < no_of_bits; i++ ) {
		*--pt_bitstream = (short)( value & MASK );
		value >>= 1;
	}
}
