#include "../include/amr_plus.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static void ExtractFormat( char *FormatString, short *ReadFlag, short *AppendFlag );
/*____________________________________________________________________________
 |
 |  FUNCTION NAME: Wave_fopen
 |____________________________________________________________________________
*/
FILE *Wave_fopen( char *Filename, char *Mode, short *NumOfChannels, long *SamplingRate, short *BitsPerSample, long *DataSize )
{
	FILE *FilePtr, *FilePtrAppend = NULL;
	short ReadFlag;
	short AppendFlag;
	char  WaveTag[4];
	long  PCM_Size = 0L;
	long  FileSize = 0L; /* Actual size will be replaces in "Wave_fclose" -procedure */
	long  FormatSize = 16L;
	short FormatTag = 1; /* Format Tag (PCM=1) */
	short BytesPerSample;
	long  BytesPerSecond;
	long  longtmp;
	short shorttmp;
	long  sRate;
	ExtractFormat( Mode, &ReadFlag, &AppendFlag );
	if( !ReadFlag ) {
		if( ( FilePtr = fopen( Filename, Mode ) ) == NULL ) {
			return FilePtr;
		}
		if( !AppendFlag ) {
			if( *BitsPerSample != 8 && *BitsPerSample != 16 ) {
				fprintf( stderr, "\n\n" );
				fprintf( stderr, "\n ERROR: Only values 16 or 8 allowed for parameter 'BitsPerSample' in function Wave_fopen" );
				fprintf( stderr, "\n        [value: %d received]", *BitsPerSample );
				fprintf( stderr, "\n\n" );
				exit( -1 );
			}
			BytesPerSample = (short)( *BitsPerSample / 8 * *NumOfChannels );
			BytesPerSecond = *SamplingRate * BytesPerSample;
			/* RIFF chunk - 12 bytes */
			fwrite( "RIFF", 1, 4, FilePtr );
			fwrite( &FileSize, 4, 1, FilePtr );
			fwrite( "WAVE", 1, 4, FilePtr );
			/* FORMAT chunk - 24 bytes */
			fwrite( "fmt ", 1, 4, FilePtr );
			fwrite( &FormatSize, 4, 1, FilePtr );
			fwrite( &FormatTag, 2, 1, FilePtr );
			fwrite( NumOfChannels, 2, 1, FilePtr );
			sRate = ( *SamplingRate );
			fwrite( &sRate, 4, 1, FilePtr );
			fwrite( &BytesPerSecond, 4, 1, FilePtr );
			fwrite( &BytesPerSample, 2, 1, FilePtr );
			fwrite( BitsPerSample, 2, 1, FilePtr );
			/* DATA chunk - Data length + 8 bytes */
			fwrite( "data", 1, 4, FilePtr );
			fwrite( &PCM_Size, 4, 1, FilePtr );
		}
		return FilePtr;
	}
	else {
		if( ( FilePtr = fopen( Filename, Mode ) ) == NULL ) {
			return FilePtr;
		}
		if( AppendFlag ) {
			FilePtrAppend = FilePtr;
		}
		fseek( FilePtr, 0, SEEK_SET );
		/* RIFF chunk - 12 bytes */
		fread( &WaveTag, 4, 1, FilePtr );
		if( strncmp( "RIFF", WaveTag, 4 ) != 0 ) {
			return NULL;
		}
		fread( &longtmp, 4, 1, FilePtr );
		fread( &WaveTag, 4, 1, FilePtr );
		if( strncmp( "WAVE", WaveTag, 4 ) != 0 ) {
			return NULL;
		}
		/* FORMAT chunk - 24 bytes */
		fread( &WaveTag, 4, 1, FilePtr );
		if( strncmp( "fmt ", WaveTag, 4 ) != 0 ) {
			return NULL;
		}
		fread( &longtmp, 4, 1, FilePtr );
		fread( &shorttmp, 2, 1, FilePtr );
		fread( NumOfChannels, 2, 1, FilePtr );
		fread( SamplingRate, 4, 1, FilePtr );
		fread( &longtmp, 4, 1, FilePtr );
		fread( &BytesPerSample, 2, 1, FilePtr );
		fread( BitsPerSample, 2, 1, FilePtr );
		/* DATA chunk - Data length + 8 bytes */
		fread( &WaveTag, 4, 1, FilePtr );
		if( strncmp( "data", WaveTag, 4 ) != 0 ) {
			return NULL;
		}
		fread( DataSize, 4, 1, FilePtr );
		*DataSize = *DataSize / BytesPerSample;
		if( AppendFlag )
			return FilePtrAppend;
		return FilePtr;
	}
}
/*____________________________________________________________________________
 |
 |  FUNCTION NAME: Wave_fclose
 |____________________________________________________________________________
*/
void Wave_fclose( FILE *FilePtr, short BitsPerSample )
{
	long DataLen;
	long FileLen;
	long StartPos;
	long EndPos;
	if( BitsPerSample != 8 && BitsPerSample != 16 ) {
		fprintf( stderr, "\n\n" );
		fprintf( stderr, "\n ERROR: Only values 16 or 8 allowed for parameter 'BitsPerSample' in function Wave_fclose" );
		fprintf( stderr, "\n\n" );
		exit( -1 );
	}
	EndPos = ftell( FilePtr );
	fseek( FilePtr, 0, SEEK_SET );
	StartPos = ftell( FilePtr );
	FileLen = ( ( EndPos - StartPos ) / ( BitsPerSample / 16 ) ) - 8;
	DataLen = FileLen - 44 + 8;
	fseek( FilePtr, 4, SEEK_SET );
	fwrite( &FileLen, 1, 4, FilePtr );
	fseek( FilePtr, 40, SEEK_SET );
	fwrite( &DataLen, 1, 4, FilePtr );
	fclose( FilePtr );
}
/*____________________________________________________________________________
 |
 |  FUNCTION NAME: ExtractFormat
 |____________________________________________________________________________
*/
static void ExtractFormat( char *FormatString, short *ReadFlag, short *AppendFlag )
{
	*ReadFlag = 0;
	*AppendFlag = 0;
	if( strchr( FormatString, (int)'a' ) != NULL )
		*AppendFlag = 1;
	if( strchr( FormatString, (int)'A' ) != NULL )
		*AppendFlag = 1;
	if( strchr( FormatString, (int)'r' ) != NULL )
		*ReadFlag = 1;
	if( strchr( FormatString, (int)'R' ) != NULL )
		*ReadFlag = 1;
}
