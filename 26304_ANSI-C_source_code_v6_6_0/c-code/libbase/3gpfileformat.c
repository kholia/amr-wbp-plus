#include "../include/amr_plus.h"
#include "ISOMovies.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VENDOR 'eric'        // 4 character vendor code "ERA/T"
#define FILETYPE_3GP6 '3gp6' // 3GPP Release 6
#define FILETYPE_3GP5                                                             \
	'3gp5' // 3GPP Release 5 (used as compatibe brand with Rel 5 when text is not \
	       // used)
#define FILETYPE_3GP4                                                             \
	'3gp4' // 3GPP Release 4 (used as compatibe brand with Rel 5 when text is not \
	       // used)

#define FILETYPE_VERSION 256 // x*256+y in Release Z.x.y of 3GPP TS 26.234

#define FRAMES_PER_SAMPLE \
	1                                 // Use -1 if we want only 1 sample containing all frames
#define MAX_FRAME_SIZE ( 2 + 4 * 80 ) // Operation code 47

// Recovered
#define MP4NewTrackIsAWB 0x10
#define MP4NewTrackIsAMRWP 0x200

// the TS ticks
static int TsTicks[13] = { 2880, 2560, 2304, 2160, 1920, 1728, 1536, 1440, 1280, 1152, 1080, 1024, 960 };

static unsigned char block_size[16] = { 18, 23, 33, 37, 41, 47, 51, 59, 61, 6, 6, 0, 0, 0, 1, 1 };
// These can be put on some state

static int TotalSampleSize = 0;
static int NumberOfFrames = 0;
static int SampleSize = 0;
static int NumberOfSamples = 0;

static MP4Handle SampleDataHandle = NULL;
static MP4Handle SampleSizeHandle = NULL;
static MP4Handle SampleDurationHandle = NULL;

static MP4Handle AMRWPEntryHandle = NULL;
static MP4Handle DecSpecInfoHandle = NULL;

static MP4Movie moov = NULL;
static MP4Track trak = NULL;
static MP4Media media = NULL;
static u32      SampleDuration;
//

// decoder reader buffer
static unsigned char *pSampleData;

static MP4TrackReader Reader;
static int            isAMRWBFileFormat;
//
MP4Err Create3GPAMRWBPlus( void )
{

	/* for the 3gp File Format*/

	MP4Err err;

	u32 TimeScale = 72000; /* which is the same as teh RTP payload format*/
	u32 objectTypeIndication = 0;
	u32 streamType = 0;
	u32 decoderBufferSize = 0;
	u32 maxBitrate = 0;
	u32 avgBitrate = 0;

	// err = MP4NewMovieWithoutIODS( &moov);
	err = ISOSetMovieBrand( moov, FILETYPE_3GP6, FILETYPE_VERSION );

	//   if (!IsExtension)
	//      err = ISOSetMovieCompatibleBrand(moov, FILETYPE_3GP5);

	/* Create new AMR track & media container, start media editing */
	err = MP4NewMovieTrack( moov, MP4NewTrackIsAMRWP, &trak );
	err = MP4NewTrackMedia( trak, &media, MP4AudioHandlerType, TimeScale, NULL );
	err = MP4BeginMediaEdits( media );

	// The following parameters are not used for AWB, set to zero.

	err = MP4SetMediaLanguage( media, "und" );

	err = MP4NewHandle( 2 * 4, &DecSpecInfoHandle );

	// setHandle32( DecSpecInfoHandle, VENDOR, 0);              //vendor
	// setHandle32( DecSpecInfoHandle, 0, 4);                   //decoder version

	err = MP4NewHandle( 0, &AMRWPEntryHandle );

	err = MP4NewSampleDescription( trak, AMRWPEntryHandle, 1, 0, 0, 0, 0, 0, DecSpecInfoHandle );

	// Create handles .....

	// Sample  Duration
	err = MP4NewHandle( 0, &SampleDurationHandle );

	SampleDuration = 0;
	// Sample Data
	err = MP4NewHandle( 0, &SampleDataHandle );
	// Sample Size
	err = MP4NewHandle( 0, &SampleSizeHandle );

	return ( err );
}

MP4Err Create3GPAMRWB( void )
{

	/* for the 3gp File Format*/

	MP4Err err;

	u32 TimeScale = 16000; /* which is the same as teh RTP payload format*/
	u32 objectTypeIndication = 0;
	u32 streamType = 0;
	u32 decoderBufferSize = 0;
	u32 maxBitrate = 0;
	u32 avgBitrate = 0;

	// err = MP4NewMovieWithoutIODS(&moov);
	err = ISOSetMovieBrand( moov, FILETYPE_3GP6, FILETYPE_VERSION );

	//   if (!IsExtension)
	err = ISOSetMovieCompatibleBrand( moov, FILETYPE_3GP5 );
	err = ISOSetMovieCompatibleBrand( moov, FILETYPE_3GP4 );
	/* Create new AMR-WB track & media container, start media editing */
	err = MP4NewMovieTrack( moov, MP4NewTrackIsAWB, &trak );
	err = MP4NewTrackMedia( trak, &media, MP4AudioHandlerType, TimeScale, NULL );
	err = MP4BeginMediaEdits( media );

	// The following parameters are not used for AWB, set to zero.

	err = MP4SetMediaLanguage( media, "und" );

	err = MP4NewHandle( 5 * 4, &DecSpecInfoHandle );

	// setHandle32( DecSpecInfoHandle, VENDOR, 0);              //vendor
	// setHandle32( DecSpecInfoHandle, 0, 4);                   //decoder version

	// setHandle32(DecSpecInfoHandle, 0x1f, 8); // mode set: default 0x1f (all
	// modes) setHandle32(DecSpecInfoHandle, 0, 12);   // mode change period
	// setHandle32(DecSpecInfoHandle, 1, 16);   // frames per sample

	err = MP4NewHandle( 0, &AMRWPEntryHandle );

	err = MP4NewSampleDescription( trak, AMRWPEntryHandle, 1, 0, 0, 0, 0, 0, DecSpecInfoHandle );

	// Create handles .....

	// Sample  Duration
	err = MP4NewHandle( 0, &SampleDurationHandle );

	SampleDuration = 0;
	// Sample Data
	err = MP4NewHandle( 0, &SampleDataHandle );
	// Sample Size
	err = MP4NewHandle( 0, &SampleSizeHandle );

	return ( err );
}

void FormatSuperFrame( EncoderConfig conf, short *Serial, int length, unsigned char *Frame, int *FrameSize )
{

	unsigned char *ptr;
	short *        Sptr;
	unsigned char  byte;
	int            i, j, k;

	ptr = Frame;
	*FrameSize = 0;
	/* write into the bitstream */
	byte = (unsigned char)conf.mode_index;
	*ptr++ = byte;
	*FrameSize = *FrameSize + 1;

	byte = (unsigned char)( conf.fscale_index & 0x1F );
	*ptr++ = byte;
	*FrameSize = *FrameSize + 1;
	for( i = 0; i < 4; i++ ) {
		Sptr = &Serial[i * ( length / 4 )];
		for( j = 0; j < length / 4 / 8; j++ ) {
			byte = 0;
			for( k = 0; k < 8; k++, Sptr++ ) {
				byte <<= 1;
				if( *Sptr != 0 )
					byte += 1;
			}

			*ptr++ = byte;
			*FrameSize = *FrameSize + 1;
		}
	}
}

void FormatAMRWBFrame( unsigned char *serial, short length, unsigned char *Frame, int *FrameSize )
{
	unsigned char *ptc_serial, *ptc_serial_old, ctemp;
	short          i;
	unsigned char  mode;
	ptc_serial = (unsigned char *)serial;
	ptc_serial_old = (unsigned char *)serial - 1;

	*Frame = ( *serial >> 3 ) << 2;
	*FrameSize = length;

	/* erease IF2 Header -> left shift 5 bits */
	for( i = 0; i < length - 1; i++ ) {
		ctemp = (unsigned char)( *ptc_serial & 0x7 ) << 5; /* save 3 lsb */
		ptc_serial++;
		ptc_serial_old++;
		ctemp += ( (unsigned char)( *ptc_serial & 0xF8 ) >> 3 ); /* save 5 msb */
		*ptc_serial_old = ctemp;
	}
	ptc_serial_old++;
	*ptc_serial_old = (unsigned char)( *ptc_serial & 0x7 ) << 5; /* save 3 lsb */

	for( i = 0; i < length; i++ ) {
		Frame[i + 1] = serial[i];
	}
}

MP4Err AddSamplesAMRWBPLUS( EncoderConfig conf, short *Serial, int length )
{
	unsigned char Frame[MAX_FRAME_SIZE];
	int           FrameSize;

	unsigned char *p;

	MP4Err        err;
	unsigned char wb_mode;

	int i;

	if( conf.extension == 1 ) { /* wb+ modes */
		FormatSuperFrame( conf, Serial, length, Frame, &FrameSize );
	}
	else {
		/* wb modes */
#ifdef IF2
		FormatAMRWBFrame( Serial, length, Frame + 1, &FrameSize );
#else
		FrameSize = length;
		p = (unsigned char *)Serial;
		for( i = 0; i < length; i++ )
			Frame[i + 1] = p[i];
#endif

		wb_mode = Frame[1] >> 3;
		Frame[0] = wb_mode;
		Frame[1] = 0;
		FrameSize++;
	}

	// Increment the SampleDataHandle Size
	err = MP4GetHandleSize( SampleDataHandle, &TotalSampleSize );
	err = MP4SetHandleSize( SampleDataHandle, TotalSampleSize + FrameSize );

	// Increment the sample duration

	if( conf.extension == 1 ) {
		SampleDuration = SampleDuration + 4 * TsTicks[conf.fscale_index - 1];
	}
	else {
		SampleDuration = SampleDuration + 1440;
	}

	// Append Data
	p = (u8 *)*SampleDataHandle;
	p += TotalSampleSize;
	for( i = 0; i < FrameSize; i++ )
		*p++ = Frame[i];

	// Increment the number of frames counter
	NumberOfFrames = NumberOfFrames + 1;
	SampleSize = SampleSize + FrameSize;
	// Check if we have enough frames stored
	if( NumberOfFrames == FRAMES_PER_SAMPLE ) {
		NumberOfSamples = NumberOfSamples + 1;
		// Resize the samples table
		err = MP4SetHandleSize( SampleSizeHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleSizeHandle, SampleSize,
		// (NumberOfSamples-1)*sizeof(u32));
		// Resize the sample duration handle
		err = MP4SetHandleSize( SampleDurationHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleDurationHandle, SampleDuration,
		// (NumberOfSamples-1)*sizeof(u32));

		SampleSize = 0;
		SampleDuration = 0;
		NumberOfFrames = 0;
	}
	return ( err );
}

MP4Err AddSamplesAMRWB( EncoderConfig conf, unsigned char *Serial, short length )
{
	unsigned char  Frame[MAX_FRAME_SIZE];
	int            FrameSize;
	int            i;
	unsigned char *p;

	MP4Err err;

	// First format the frame, needed only if the IF2 format is selected
#ifdef IF2
	FormatAMRWBFrame( Serial, length, Frame, &FrameSize );
#else
	FrameSize = length;

	for( i = 0; i < length; i++ )
		Frame[i] = Serial[i];
#endif

	// Frame = Serial;
	// FrameSize = length;
	// Increment the SampleDataHandle Size
	err = MP4GetHandleSize( SampleDataHandle, &TotalSampleSize );
	err = MP4SetHandleSize( SampleDataHandle, TotalSampleSize + FrameSize );

	// Increment the sample duration

	SampleDuration = SampleDuration + 320;

	// Append Data
	p = (u8 *)*SampleDataHandle;
	p += TotalSampleSize;
	for( i = 0; i < FrameSize; i++ )
		*p++ = Frame[i];

	// Increment the number of frames counter
	NumberOfFrames = NumberOfFrames + 1;
	SampleSize = SampleSize + FrameSize;
	// Check if we have enough frames stored
	if( NumberOfFrames == FRAMES_PER_SAMPLE ) {
		NumberOfSamples = NumberOfSamples + 1;
		// Resize the samples table
		err = MP4SetHandleSize( SampleSizeHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleSizeHandle, SampleSize,
		// (NumberOfSamples-1)*sizeof(u32)); Resize the sample duration handle
		err = MP4SetHandleSize( SampleDurationHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleDurationHandle, SampleDuration,
		// (NumberOfSamples-1)*sizeof(u32));

		SampleSize = 0;
		SampleDuration = 0;
		NumberOfFrames = 0;
	}
	return ( err );
}

MP4Err WriteSamplesAMRWBPlus( EncoderConfig conf, void *Serial, int length )
{

	/* test if backward compatible file format */
	if( conf.bc == 1 ) {

		AddSamplesAMRWB( conf, (unsigned char *)Serial, (short)length );
	}
	else {

		AddSamplesAMRWBPLUS( conf, (short *)Serial, (int)length );
	}
	return 0;
}
MP4Err Close3GP( char *filename )
{

	MP4Err err;
	u64    Duration;
	// This is the end, resize the last sample

	if( NumberOfFrames != 0 ) {
		NumberOfSamples = NumberOfSamples + 1;
		// Resize
		err = MP4SetHandleSize( SampleSizeHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleSizeHandle, SampleSize,
		// (NumberOfSamples-1)*sizeof(u32));

		// Resize the sample duration handle
		err = MP4SetHandleSize( SampleDurationHandle, NumberOfSamples * sizeof( u32 ) );
		// setHandle32(SampleDurationHandle, SampleDuration,
		// (NumberOfSamples-1)*sizeof(u32));

		SampleSize = 0;
		SampleDuration = 0;
		NumberOfFrames = 0;
	}

	err = MP4AddMediaSamples( media, SampleDataHandle, NumberOfSamples, SampleDurationHandle, SampleSizeHandle, AMRWPEntryHandle, NULL, NULL );

	err = MP4EndMediaEdits( media );
	err = MP4GetMediaDuration( media, &Duration );
	err = MP4InsertMediaIntoTrack( trak, 0, 0, Duration, 1 );

	err = MP4WriteMovieToFile( moov, filename );
	err = MP4DisposeMovie( moov );

	// Dispose of all handles
	err = MP4DisposeHandle( SampleSizeHandle );
	err = MP4DisposeHandle( SampleDataHandle );
	err = MP4DisposeHandle( SampleDurationHandle );
	err = MP4DisposeHandle( SampleSizeHandle );
	err = MP4DisposeHandle( AMRWPEntryHandle );
	return ( err );
}

MP4Err GetBrand( MP4Movie moov, int *isom_comp, int *rel4_comp, int *rel5_comp )
{
	MP4Err err;
	u32    brand, minorversion;
	u8     s1, s2, s3, s4;

	err = ISOGetMovieBrand( moov, &brand, &minorversion );
	if( err ) {
		fprintf( stderr, "\nMovie brand not specified\n" );
		return err;
	}

	s1 = ( brand >> 24 ) & 0xff;
	s2 = ( brand >> 16 ) & 0xff;
	s3 = ( brand >> 8 ) & 0xff;
	s4 = brand & 0xff;

	printf( "\nFile type brand '%c%c%c%c' minor version %i\n", s1, s2, s3, s4, minorversion );

	*rel4_comp = ( ISOIsMovieCompatibleBrand( moov, '3gp4' ) != 0 );
	*rel5_comp = ( ISOIsMovieCompatibleBrand( moov, '3gp5' ) != 0 );
	*isom_comp = ( ISOIsMovieCompatibleBrand( moov, 'isom' ) != 0 );

	printf( "Compatible brands: 'isom' (" );
	if( *isom_comp ) {
		printf( "yes" );
	}
	else {
		printf( "no" );
	}

	printf( ") '3gp4' (" );
	if( *rel4_comp ) {
		printf( "yes" );
	}
	else {
		printf( "no" );
	}
	printf( ") '3gp5' (" );

	if( *rel5_comp ) {
		printf( "yes" );
	}
	else {
		printf( "no" );
	}
	printf( ")\n" );

	return err;
}

MP4Err FindWBPlusTrack( MP4Movie moov, int TrackNumber, u32 *trackType )
{
	MP4Err    err;
	MP4Track  track;
	MP4Media  media;
	MP4Handle sampleDescriptionH;
	u32       handlerType;
	u32       sampleEntryAtomType;
	u32       timeScale;
	u32       trackID;
	u64       duration;
	u8        s1, s2, s3, s4;
	u8 *      q;

	err = MP4GetMovieIndTrack( moov, TrackNumber, &track );
	if( err )
		return err;

	err = MP4GetTrackMedia( track, &media );
	if( err )
		return err;

	// get track ID
	err = MP4GetTrackID( track, &trackID );
	if( err )
		return err;

	// get time scale and duration
	err = MP4GetMediaTimeScale( media, &timeScale );
	if( err )
		return err;

	err = MP4GetMediaDuration( media, &duration );
	if( err )
		return err;

	// get handler type
	err = MP4GetMediaHandlerDescription( media, &handlerType, NULL );
	if( err )
		return err;

	// get sample description
	err = MP4NewHandle( 0, &sampleDescriptionH );
	if( err )
		return err;

	err = MP4GetMediaSampleDescription( media, 1, sampleDescriptionH, NULL );
	if( err )
		return err;

	q = (u8 *)( 4 + *sampleDescriptionH );
	s1 = *q++;
	s2 = *q++;
	s3 = *q++;
	s4 = *q++;
	sampleEntryAtomType = ( s1 << 24 ) + ( s2 << 16 ) + ( s3 << 8 ) + s4;

	printf( "Track %i '%c%c%c%c' ", trackID, s1, s2, s3, s4 );

	err = MP4NotImplementedErr;

	if( handlerType == MP4AudioHandlerType ) {
		/* found an audio trak */
		if( sampleEntryAtomType == 'sawb' ) {

			*trackType = MP4NewTrackIsAWB;
			err = MP4NoErr;
		}
		else if( sampleEntryAtomType == 'sawp' ) {
			*trackType = MP4NewTrackIsAMRWP;
			err = MP4NoErr;
		}
	}
	printf( " Duration %.3fs  Time scale %d\n", (double)duration / timeScale, timeScale );

	MP4DisposeHandle( sampleDescriptionH );
	return err;
}

// MP4Err GetNextFrame3GP(short *extension, short *mode, short *st_mode, short *fst, void *serial, int init) {
MP4Err GetNextFrame3GP( short *tfi, int *bfi, short *extension, short *mode, short *st_mode, short *fst, void *serial, int init )
{
	MP4Err err;
	u32    outSampleFlags;
	u32    outCTS;
	u32    outDTS;
	short *ptr;
	int    frame_size;
	int    mode_index, j, k;
	int    fscale_index, index;
	// int tfi, fst_index, nb_bits;
	int            fst_index, nb_bits;
	unsigned char  byte;
	int            sample_size;
	unsigned char  ctemp;
	int            nb_bytes;
	int            i;
	unsigned char *ptc_serial, *ptc_serial_new;
	short *        pts_serial;

	if( SampleSize == 0 ) {
		// read in a new sample
		err = MP4TrackReaderGetNextAccessUnit( Reader, SampleDataHandle, &SampleSize, &outSampleFlags, &outCTS, &outDTS );
		if( err ) {
			// Dispose of the handles and exit
			ISODisposeTrackReader( Reader );
			MP4DisposeHandle( SampleDataHandle );
			return ( err );
		}

		pSampleData = ( *SampleDataHandle );
	}

	sample_size = SampleSize;

	// at this point sample size should contain something

	if( isAMRWBFileFormat ) {
		mode_index = *pSampleData;
		mode_index = 15 & ( mode_index >> 3 );
		*extension = 0;
		*st_mode = -1;
		*fst = 0;
		*mode = mode_index;

		if( !init ) {

			frame_size = block_size[*mode];

			ptc_serial = ( (unsigned char *)serial ); // leave 1 byte for mms header
			for( j = 0; j < frame_size; j++ ) {
				*ptc_serial++ = *pSampleData++;
				SampleSize--;
			}
			//

			// add IF2 header
		}
	}
	else {
		// read mode index
		mode_index = *pSampleData++;
		SampleSize--;

		// read isf index
		fscale_index = *pSampleData++;
		SampleSize--;

		// get the mode and isf indices

		if( mode_index > 47 || mode_index < 0 ||     /* mode unknown */
		    mode_index == 14 || mode_index == 15 ||  /* Frame lost or ereased */
		    ( mode_index == 9 && *extension == 1 ) ) /* WB SID in WB+ frame */
		{
			*tfi = 0; // note

			/* If mode_index unknown, frame is ereased  or NO_DATA */
			if( mode_index == 14 ) {
				*mode = 14;
				bfi[0] = 1; // note
				bfi[1] = 1;
				bfi[2] = 1;
				bfi[3] = 1;
			}
			else {
				*mode = 15;
				bfi[0] = 0; // note
				bfi[1] = 0;
				bfi[2] = 0;
				bfi[3] = 0;
			}
		}
		else if( mode_index > 15 ) /* wb+ */
		{
			if( mode_index < 24 ) /* Mono mode only */
			{
				*mode = mode_index - 16;
				*st_mode = -1;
			}
			else {
				index = mode_index - 24;
				*mode = miMode[2 * index];
				*st_mode = miMode[2 * index + 1];
			}
			*extension = 1;

			// tfi = (fscale_index & 0xc0);
			// *tfi = (fscale_index & 0xc0);
			*tfi = 0; // note
			bfi[0] = 0;
			bfi[1] = 0;
			bfi[2] = 0;
			bfi[3] = 0;

			fst_index = ( fscale_index & 0x1F );
			if( fst_index < 1 )
				fst_index = 1; /* prevent isf < 0.5 */
			*fst = isfIndex[fst_index];
		}
		else { /* wb in wb+ file format */
			/* special AMR-WB+ modes */
			if( mode_index == 10 ) {
				*extension = 1;
				*mode = 2; /* 14m */
				*st_mode = -1;
			}
			else if( mode_index == 11 ) {
				*extension = 1;
				*mode = 2; /* 18s */
				*st_mode = 6;
			}
			else if( mode_index == 12 ) {
				*extension = 1;
				*mode = 7; /* 24m */
				*st_mode = -1;
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

			// tfi = (fscale_index & 0xc0);
			// *tfi = (fscale_index & 0xc0);
			*tfi = 0; // note
			bfi[0] = 0;
			bfi[1] = 0;
			bfi[2] = 0;
			bfi[3] = 0;
			fst_index = ( fscale_index & 0x1F );
			if( fst_index != 0 && fst_index != 8 ) {
				fprintf( stderr,
				    "Internal Sampling Frequency not supported with AMW WB "
				    "and caracterized WB+ modes " );
				exit( EXIT_FAILURE );
			}
			*fst = isfIndex[fst_index];
		}

		if( init ) {
			SampleSize += 2;
			pSampleData -= 2;
		}
		else {
			if( *extension == 1 ) {
				// get the frame size in bytes
				nb_bits = NBITS_CORE[*mode] + NBITS_BWE;
				if( *st_mode >= 0 ) {
					nb_bits += ( StereoNbits[*st_mode] + NBITS_BWE );
				}
				frame_size = nb_bits / 8;

				ptr = (short *)serial;
				for( j = 0; j < frame_size; j++ ) {
					byte = *pSampleData++;
					SampleSize--;
					for( k = 0; k < 8; k++, ptr++ ) {
						*ptr = ( byte & (short)128 ) == (short)128;
						byte <<= 1;
					}
				}
			}
			else {
#ifdef IF2
				/* not implemented */
#else
				// nb_bits= NBITS_CORE_AMR_WB[*mode]/4;
				frame_size = block_size[*mode];

				ptc_serial = ( (unsigned char *)serial ) + 1; // leave 1 byte for mms header
				for( j = 0; j < frame_size - 1; j++ ) {
					*ptc_serial++ = *pSampleData++;
					SampleSize--;
				}

				ptc_serial = ( (unsigned char *)serial );
				*ptc_serial = ( *mode << 3 ) + 4;

#endif
			}
		}
	}
	return MP4NoErr;
}

/*
   Open a 3GP file for reading
*/
MP4Err Open3GP_( char *filename, int verbose, DecoderConfig *conf )
{

	/* for the 3gp File Format*/

	MP4Err err;
	int    OpenMovieFlags;
	int    TrackCount;
	int    isom_comp;
	int    rel4_comp;
	int    rel5_comp;
	int    i;
	u32    trackType;
	int    AWBFound = 0, AMRWPFound = 0;
	int    AWBtrack, AMRWPtrack;
	int    trackNumber;
	// open movie (3GP file)
	if( verbose ) {
		OpenMovieFlags = MP4OpenMovieDebug;
	}
	else {
		OpenMovieFlags = MP4OpenMovieNormal;
	}

	err = MP4OpenMovieFile( &moov, filename, OpenMovieFlags );
	if( err ) {
		return ( err );
	}

	err = GetBrand( moov, &isom_comp, &rel4_comp, &rel5_comp );

	// go through tracks of movie
	err = MP4GetMovieTrackCount( moov, &TrackCount );

	for( i = 1; i <= (int)TrackCount; i++ ) {
		err = FindWBPlusTrack( moov, i, &trackType );
		if( err ) {
			fprintf( stderr, "\n Unable de decode" );
			return ( err );
		}
		if( ( trackType & MP4NewTrackIsAWB ) ) {
			if( !AWBFound ) {
				AWBFound = 1;
				AWBtrack = i;
			}
			else
				printf( "More than one AWB track found. This track will be ignored!\n" );
		}
		if( ( trackType & MP4NewTrackIsAMRWP ) ) {
			if( !AMRWPFound ) {
				AMRWPFound = 1;
				AMRWPtrack = i;
			}
			else
				printf(
				    "More than one AWB-WB+ track found. This track will be ignored!\n" );
		}
	}

	if( AWBFound && AMRWPFound ) {
		printf( "\n AMR-WB and AMR-WB+ tracks found. Decoding only WB+ tracks" );
		AWBFound = 1;
	}

	if( AWBFound )
		trackNumber = AWBtrack;
	if( AMRWPFound )
		trackNumber = AMRWPtrack;

	err = MP4GetMovieIndTrack( moov, trackNumber, &trak );
	if( err ) {
		return ( err );
	}
	err = MP4GetTrackMedia( trak, &media );
	if( err ) {
		return ( err );
	}

	// err = MP4GetMediaHandlerDescription( media, &HandlerType, NULL );
	if( err ) {
		return ( err );
	}
	err = MP4CreateTrackReader( trak, &Reader );
	if( err ) {
		return ( err );
	}

	// err = reportAMRSpecificInfo(trak, WB);
	err = MP4NewHandle( 0, &SampleDataHandle );

	if( AMRWPFound ) {
		conf->extension = 1;
		isAMRWBFileFormat = 0;
	}
	else {
		conf->extension = 0;
		isAMRWBFileFormat = 1;
	}

	/* GetNextFrame3GP(&(conf->extension), &(conf->mode), &(conf->st_mode),
                  &(conf->fscale), NULL, 1); */
	return err;
}

/*
   Open a 3GP file for reading
*/
MP4Err Open3GP( short *tfi, int *bfi, char *filename, int verbose, DecoderConfig *conf )
{

	/* for the 3gp File Format*/

	MP4Err err;
	int    OpenMovieFlags;
	int    TrackCount;
	int    isom_comp;
	int    rel4_comp;
	int    rel5_comp;
	int    i;
	u32    trackType;
	int    AWBFound = 0, AMRWPFound = 0;
	int    AWBtrack, AMRWPtrack;
	int    trackNumber;
	// open movie (3GP file)
	if( verbose ) {
		OpenMovieFlags = MP4OpenMovieDebug;
	}
	else {
		OpenMovieFlags = MP4OpenMovieNormal;
	}

	err = MP4OpenMovieFile( &moov, filename, OpenMovieFlags );
	if( err ) {
		return ( err );
	}

	err = GetBrand( moov, &isom_comp, &rel4_comp, &rel5_comp );

	// go through tracks of movie
	err = MP4GetMovieTrackCount( moov, &TrackCount );

	for( i = 1; i <= (int)TrackCount; i++ ) {
		err = FindWBPlusTrack( moov, i, &trackType );
		if( err ) {
			fprintf( stderr, "\n Unable de decode" );
			return ( err );
		}
		if( ( trackType & MP4NewTrackIsAWB ) ) {
			if( !AWBFound ) {
				AWBFound = 1;
				AWBtrack = i;
			}
			else
				printf( "More than one AWB track found. This track will be ignored!\n" );
		}
		if( ( trackType & MP4NewTrackIsAMRWP ) ) {
			if( !AMRWPFound ) {
				AMRWPFound = 1;
				AMRWPtrack = i;
			}
			else
				printf(
				    "More than one AWB-WB+ track found. This track will be ignored!\n" );
		}
	}

	if( AWBFound && AMRWPFound ) {
		printf( "\n AMR-WB and AMR-WB+ tracks found. Decoding only WB+ tracks" );
		AWBFound = 1;
	}

	if( AWBFound )
		trackNumber = AWBtrack;
	if( AMRWPFound )
		trackNumber = AMRWPtrack;

	err = MP4GetMovieIndTrack( moov, trackNumber, &trak );
	if( err ) {
		return ( err );
	}
	err = MP4GetTrackMedia( trak, &media );
	if( err ) {
		return ( err );
	}

	// err = MP4GetMediaHandlerDescription( media, &HandlerType, NULL );
	if( err ) {
		return ( err );
	}
	err = MP4CreateTrackReader( trak, &Reader );
	if( err ) {
		return ( err );
	}

	// err = reportAMRSpecificInfo(trak, WB);
	err = MP4NewHandle( 0, &SampleDataHandle );

	if( AMRWPFound ) {
		conf->extension = 1;
		isAMRWBFileFormat = 0;
	}
	else {
		conf->extension = 0;
		isAMRWBFileFormat = 1;
	}

	GetNextFrame3GP( tfi, bfi, &( conf->extension ), &( conf->mode ), &( conf->st_mode ), &( conf->fscale ), NULL, 1 );

	/* GetNextFrame3GP(&(conf->extension), &(conf->mode), &(conf->st_mode),
                  &(conf->fscale), NULL, 1); */
	return err;
}
