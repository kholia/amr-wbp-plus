#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

extern const UWord8 block_size[];
#define NO_DATA -3
/*-------------------------------------------------------------------------------------*
 * Funtion  WriteHeader()
 * ~~~~~~~~~~~~~~~~~~
 * Write AMR-WB and AMR-WB+ header according to payload selected (see 3gpp ts26.304)
 *
 *                            AMR-WB+
 *
 * Raw File Format (FRTP)
 *  | M.Index (7 bits) | T.F. Index (2 bits) | I.S.F Index (5 bits) | AMR-WB+ bitstream...
 *
 *                            AMR-WB
 *
 *  | Mono Rate (8 bits) | .........AMR-WB bitstream in IF2 Format (without IF2 Header)...
 *--------------------------------------------------------------------------------------*/

void WriteHeader( EncoderConfig conf, short length, short offset, FILE *f_serial )
{

	short         codec_mode;
	unsigned char byte;

	if( conf.extension > 0 ) {
		/* write into the bitstream */
		byte = (unsigned char)conf.mode_index;
		fwrite( &byte, sizeof( char ), 1, f_serial ); /*write mode index information */

		byte = ( offset << 6 );
		byte += (unsigned char)( conf.fscale_index & 0x1F );
		fwrite( &byte, sizeof( char ), 1, f_serial );
	}
	else /* AMR WB */
	{
		codec_mode = conf.mode;
		if( length == 1 )
			codec_mode = 15;
		else if( length == 6 )
			codec_mode = 9;
		byte = (unsigned char)codec_mode;
		fwrite( &byte, sizeof( char ), 1, f_serial );

		byte = 0;
		byte = ( offset << 6 );
		fwrite( &byte, sizeof( char ), 1, f_serial );
	}
}
/*-----------------------------------------------------------------*
 * Funtion  WriteBitstreamPlus()
 * ~~~~~~~~~~~~~~~~~~
 * Write AMR-WB+ bitstream
 *-----------------------------------------------------------------*/

void WriteBitstreamPlus( EncoderConfig conf, short length, short offset, short *serial, FILE *f_serial )
{
	unsigned char byte;
	short         j, k, nb_byte, *pt_serial, *ptr;

	pt_serial = (short *)serial;
	nb_byte = ( ( length / 4 ) + 7 ) / 8;
	ptr = &pt_serial[offset * ( length / 4 )];
	for( j = 0; j < nb_byte; j++ ) {
		byte = 0;
		for( k = 0; k < 8; k++, ptr++ ) {
			byte <<= 1;
			if( *ptr != 0 )
				byte += 1;
		}
		fwrite( &byte, sizeof( char ), 1, f_serial );
	}
}
/*-----------------------------------------------------------------*
 * Funtion  WriteBitstream()
 * ~~~~~~~~~~~~~~~~~~
 * Write AMR-WB bitstream
 *-----------------------------------------------------------------*/

void WriteBitstream( EncoderConfig conf, short length, short offset, unsigned char *serial, FILE *f_serial )
{
	unsigned char *ptc_serial, *ptc_serial_old;
	short          i;

	/* Remove IF2 header */

	ptc_serial = (unsigned char *)serial + 1;
	ptc_serial_old = (unsigned char *)serial;
	for( i = 0; i < length - 1; i++ ) {
		*ptc_serial_old = (unsigned char)( *ptc_serial );
		ptc_serial++;
		ptc_serial_old++;
	}
	fwrite( serial, sizeof( unsigned char ), length - 1, f_serial );
}
/*-------------------------------------------------------------------------------------*
 * Funtion  ReadHeader()
 * ~~~~~~~~~~~~~~~~~~
 * Read AMR-WB and AMR-WB+ header according to payload selected (see 3gpp ts26.304)
 *
 *                            AMR-WB+
 * Raw File Format (FRAW)
 *  | M.Index (7 bits) | T.F. Index (2 bits) | I.S.F Index (5 bits) | AMR-WB+ bitstream...
 *
 *                            AMR-WB
 *
 *  | Mono Rate (8 bits) | .........AMR-WB bitstream in IF2 Format (without IF2 Header)...
 *--------------------------------------------------------------------------------------*/

short ReadHeader( short *tfi, int *bfi, short FileFormat, short *extension, short *mode, short *st_mode, short *fst, short offset, FILE *f_serial )
{
	short         mode_index, fst_index, index, nb_read = 0;
	unsigned char byte;

	nb_read += fread( &byte, sizeof( unsigned char ), 1, f_serial );
	if( nb_read == 0 )
		return nb_read;

	mode_index = ( byte & 127 );

	/* If frame ereased : don't change old conf just modify mode */
	if( mode_index > 47 || mode_index < 0 ||     /* mode unknown */
	    mode_index == 14 || mode_index == 15 ||  /* Frame lost or ereased */
	    ( mode_index == 9 && *extension == 1 ) ) /* WB SID in WB+ frame  not supported case so declare a NO_DATA*/
	{

		nb_read += fread( &byte, sizeof( unsigned char ), 1, f_serial ); /* read one more byte to ensure empty header */
		*tfi = ( byte & 0xc0 ) >> 6;                                     /* tfi extrapolated by RTP packetizer */
		fst_index = ( byte & 0x1F );
		*fst = isfIndex[fst_index];

		if( mode_index == 14 ) /* frame lost WB or WB+*/
		{
			bfi[*tfi] = 1;
		}
		else if( mode_index == 15 ) {
			bfi[*tfi] = 0; /* DTX in WB reset BFI vector */
		}
		*mode = mode_index;
		return NO_DATA; /* There is no more data to read */
	}
	*st_mode = -1;

	if( mode_index > 15 ) /* wb+ */
	{
		if( mode_index < 24 ) /* Mono mode only */
		{
			*mode = mode_index - 16;
		}
		else {
			index = mode_index - 24;
			*mode = miMode[2 * index];
			*st_mode = miMode[2 * index + 1];
		}
		*extension = 1;
		nb_read += fread( &byte, sizeof( unsigned char ), 1, f_serial );
		*tfi = ( byte & 0xc0 ) >> 6;
		fst_index = ( byte & 0x1F );
		if( fst_index < 1 )
			fst_index = 1; /* prevent isf < 0.5 */
		*fst = isfIndex[fst_index];
	}
	else /* WB and caracterize WB+*/
	{
		if( mode_index == 10 ) {
			*extension = 1;
			*mode = 2; /* 14m */
		}
		else if( mode_index == 11 ) {
			*extension = 1;
			*mode = 2; /* 18s */
			*st_mode = 6;
		}
		else if( mode_index == 12 ) {
			*extension = 1;
			*mode = 7; /* 24m */
		}
		else if( mode_index == 13 ) {
			*extension = 1;
			*mode = 5; /* 24s */
			*st_mode = 7;
		}

		else {
			*extension = 0;
			*mode = mode_index;
		}
		nb_read += fread( &byte, sizeof( unsigned char ), 1, f_serial );
		*tfi = ( byte & 0xc0 ) >> 6;
		fst_index = ( byte & 0x1F );
		if( fst_index != 0 && fst_index != 8 ) {
			fprintf( stderr, "Internal Sampling Frequency not supported with AMW WB and caracterized WB+ modes " );
			exit( EXIT_FAILURE );
		}
		*fst = isfIndex[fst_index];
	}
	bfi[*tfi] = 0; /* Good frame */

	return nb_read;
}
/*-----------------------------------------------------------------*
 * Funtion  WriteBitstreamPlus()
 * ~~~~~~~~~~~~~~~~~~
 * Write AMR-WB+ bitstream
 *-----------------------------------------------------------------*/

static short ReadBitstreamPlus( short nb_bits, short nb_byte, short *serial, FILE *f_serial, short offset )
{
	unsigned char byte;
	short         j, k, n, *ptr;

	ptr = &serial[offset * ( nb_bits / 4 )];
	n = 0;
	for( j = 0; j < nb_byte; j++ ) {
		n += fread( &byte, sizeof( unsigned char ), 1, f_serial );
		for( k = 0; k < 8; k++, ptr++ ) {
			*ptr = ( byte & (short)128 ) == (short)128;
			byte <<= 1;
		}
	}

	return n;
}
/*-----------------------------------------------------------------*
 * Funtion  WriteBitstream()
 * ~~~~~~~~~~~~~~~~~~
 * Write AMR-WB bitstream + creation of IF2 Header
 *-----------------------------------------------------------------*/

static short ReadBitstream( short nb_byte, unsigned char *serialAmrwb, short mode, FILE *f_serial )
{
	unsigned char *ptc_serial, ctemp;
	unsigned char *ptc_serial_new;

	short n, i;
	/* update mode changes */
	n = fread( serialAmrwb, sizeof( unsigned char ), nb_byte - 1, f_serial );

	ptc_serial = serialAmrwb + nb_byte - 2;
	ptc_serial_new = serialAmrwb + nb_byte - 1;

	for( i = 0; i < nb_byte - 1; i++ ) {
		*ptc_serial_new = *ptc_serial;

		ptc_serial--;
		ptc_serial_new--;
	}
	/* add IF2 Header */
	ctemp = (unsigned char)( 1 << 2 );     /* Add FQI        */
	ctemp += (unsigned char)( mode << 3 ); /* Add Frame Type */

	*ptc_serial_new = ctemp;

	return n;
}

short ReadRawFile( short *tfi, int *bfi, DecoderConfig *conf, short *extension, short *mode, short *st_mode, short *fst, FILE *f_serial, void *serial )
{
	short i, n = 0, nb_bits, nb_byte, *pt_serial, *ptr;
	pt_serial = (short *)serial;
	for( i = 0; i < 4; i++ ) {
		n = ReadHeader( tfi, bfi, conf->FileFormat, extension, mode, st_mode, fst, i, f_serial );
		if( !n )
			break;

		/* assume there is no amrwb -> wb+ switching inside superframe */
		if( n != NO_DATA ) {
			/* update mode and st_mode only if mode_index != (14 ||15) */
			conf->mode = *mode;
			conf->st_mode = *st_mode;
			if( *extension > 0 ) {
				nb_bits = get_nb_bits( *extension, *mode, *st_mode );
				ptr = &pt_serial[i * ( nb_bits / 4 )];

				nb_byte = ( ( nb_bits / 4 ) + 7 ) / 8;
				n = ReadBitstreamPlus( nb_bits, nb_byte, pt_serial, f_serial, i );
				if( n != nb_byte )
					break;
			}
			else {
				nb_byte = block_size[*mode];
				n = ReadBitstream( nb_byte, (unsigned char *)serial, *mode, f_serial );
				break;
			}
		}
		else if( *extension == 0 ) {
			if( *mode == 15 ) {                        /* DTX FRAME (NO DATA)*/
				( (unsigned char *)serial )[0] = 0x7C; /* need in AMR WB */
			}
			else if( *mode == 14 )                     /* Frame lost */
				( (unsigned char *)serial )[0] = 0x74; /* need in AMR WB */
			break;
		}
	}
	return n;
}
