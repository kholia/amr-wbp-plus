#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/amr_plus.h"
#include "../include/wbplus3gplib.h"
#include "../lib_amr/dec_if.h"
#include "../lib_amr/enc_if.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

extern const UWord8 block_size[];

static void copyright( void )
{
	fprintf( stderr, "\n" );
}
static void usage( char *argv )
{
	fprintf( stderr,
	    "Usage: %s [-ff <3gp/raw> -fs <fs_Hz>][-mono] [-limiter] -if "
	    "<infile.wb+> -of <outfile.wav> -fer <error.txt>\n",
	    argv );
	fprintf( stderr, "\n" );
	fprintf( stderr, "-ff     file format used (default is 3gp)\n" );
	fprintf( stderr, "-fs     frequency sampling\n" );
	fprintf( stderr, "-mono   mono output \n" );
	fprintf( stderr, "-limiter    avoids output clipping (recommended)\n" );
	fprintf( stderr, "-if     input AMRWB+ bitstream file\n" );
	fprintf( stderr, "-of     output audio WAV file\n" );
	fprintf( stderr, "-fer    frame erasures file (0/1)" );
	fprintf( stderr, "\n" );
}

static float get_bitrate( DecoderConfig *conf )
{
	if( conf->fscale != 0 ) {
		return (float)( get_nb_bits( conf->extension, conf->mode, conf->st_mode ) * conf->fscale ) / ( 80.0f * FSCALE_DENOM );
	}
	else {
		return (float)get_nb_bits( conf->extension, conf->mode, conf->st_mode ) / 80.0f;
	}
}
static void parsecmdline( int argc, char *argv[], char **input_filename, char **output_filename, char **fer_filename, DecoderConfig *conf )
{
	conf->fs = 0;
	conf->mono_dec_stereo = 0;
	conf->limiter_on = 0;
	conf->fer_sim = 0;
	conf->FileFormat = F3GP;
	if( argc < 5 ) {
		usage( *argv );
		exit( EXIT_FAILURE );
	}
	argc--;
	argv++;
	while( argc > 0 ) {
		if( !strcmp( *argv, "-fs" ) ) {
			argv++;
			argc--;
			conf->fs = atoi( *argv );
		}
		else if( !strcmp( *argv, "-ff" ) ) {
			argv++;
			argc--;
			if( !strcmp( *argv, "3gp" ) ) {
				conf->FileFormat = F3GP;
			}
			else {
				conf->FileFormat = FRAW;
			}
		}

		else if( !strcmp( *argv, "-mono" ) ) {
			conf->mono_dec_stereo = 1;
		}
		else if( !strcmp( *argv, "-limiter" ) ) {
			conf->limiter_on = 1;
		}
		else if( !strcmp( *argv, "-if" ) ) {
			argv++;
			argc--;
			*input_filename = *argv;
		}
		else if( !strcmp( *argv, "-of" ) ) {
			argv++;
			argc--;
			*output_filename = *argv;
		}
		else if( !strcmp( *argv, "-fer" ) ) {
			argv++;
			argc--;
			*fer_filename = *argv;
			conf->fer_sim = 1;
		}
		else {
			fprintf( stderr, "Unknown option %s\n", *argv );
			exit( EXIT_FAILURE );
		}
		argv++;
		argc--;
	}
}

static void select_fs( DecoderConfig *conf )
{
	/* default sampling rate if undefined */
	if( conf->fs == 0 ) {
		if( conf->extension == 0 )
			conf->fs = 16000;
		else if( conf->fscale == 0 )
			conf->fs = 24000;
		else
			conf->fs = 48000;
	}
	/* user specified sampling rate */
	else {
		if( conf->extension == 0 ) {
			if( conf->fs != 16000 ) {
				fprintf( stderr, "Sampling rate not supported" );
				exit( EXIT_FAILURE );
			}
		}
		else {
			if( conf->fscale == 0 ) {
				if( ( conf->fs != 16000 ) && ( conf->fs != 24000 ) && ( conf->fs != 8000 ) ) {
					fprintf( stderr, "Sampling rate not supported" );
					exit( EXIT_FAILURE );
				}
			}
			else {
				if( ( conf->fs != 44100 ) && ( conf->fs != 48000 ) ) {
					fprintf( stderr, "Sampling rate not supported" );
					exit( EXIT_FAILURE );
				}
			}
		}
	}
}
static void set_frame_length( DecoderConfig *conf, int *L_frame )
{
	switch( conf->fs ) {
	case 8000:
		*L_frame = L_FRAME8k;
		break;
	case 16000:
		*L_frame = L_FRAME16kPLUS;
		break;
	case 24000:
		*L_frame = L_FRAME24k;
		break;
	case 32000:
		*L_frame = L_FRAME32k;
		break;
#ifdef FILTER_44kHz
	case 44100:
		*L_frame = L_FRAME44k;
		break;
#endif
#ifdef FILTER_48kHz
	case 48000:
		*L_frame = L_FRAME48k;
		break;
#endif
	default:
		fprintf( stderr, "Sampling rate not supported" );
		exit( 1 );
	}
}
static void interleave( float *right, float *left, float *out, int length )
{
	int i;

	for( i = length - 1; i >= 0; i-- ) {
		out[( i * 2 ) + 1] = right[i];
		out[i * 2] = left[i];
	}
}
static void simple_frame_limiter( float *signal, float *gain, int lg )
{
	int   i;
	int   pos_start, pos_end;
	float amp, slope;
	float max_amp;
	float target_gain;

	max_amp = -1;
	pos_start = -1;
	for( i = 0; i < lg; i++ ) {
		amp = (float)fabs( signal[i] );
		if( amp > 32767.0f ) {
			if( pos_start < 0 )
				pos_start = i;
			pos_end = i;
			if( amp > max_amp ) {
				max_amp = amp;
			}
		}
	}
	if( max_amp > 0 ) {
		target_gain = 32767.0f / max_amp;
		slope = ( target_gain - *gain ) / ( (float)pos_start + 1.0f );
		for( i = 0; i < lg; i++ ) {
			if( i < pos_start + 1 )
				*gain += slope;
			signal[i] *= *gain;
		}
	}
	else {
		if( *gain != 1.0 ) {
			slope = ( 1.0f - *gain ) / ( (float)lg );
			for( i = 0; i < lg; i++ ) {
				*gain += slope;
				signal[i] *= *gain;
			}
		}
	}
}

void main( int argc, char *argv[] )
{
	FILE *f_serial; /* File of serial bits for transmission  */
	FILE *f_syn;    /* File of synthesis speech              */
	FILE *f_fer;    /* File of frame errasures              */

	short               speech16[L_FRAME48k * 2];
	Decoder_State_Plus *st_d;
	void *              stAmrwbDec;
	unsigned short      serial[NBITS_MAX]; /* serial parameters.                  */

	short         bitsPerSample = 16;
	int           bfi[4] = { 0, 0, 0, 0 };
	int           L_frame;
	long          frame, dataSize;
	float         channel_right[4 * L_FRAME48k];
	float         channel_left[2 * L_FRAME48k];
	float         mem_down_right[2 * L_FILT_DECIM_FS], mem_down_left[2 * L_FILT_DECIM_FS];
	int           frac_down_right, frac_down_left;
	int           fac_up, fac_down, nb_samp_fs;
	long          fs_output;
	short         nb_samp;
	short         mode, extension;
	short         st_mode;
	short         num_channels;
	char *        input_filename = "";
	char *        output_filename = "";
	char *        fer_filename = "";
	short         fst_switched;
	short         old_st_mode;
	short         fst;
	DecoderConfig conf;
	float         left_limiter_gain;
	float         right_limiter_gain;
	short         i, j;

	short tfi = 0;
	float old_bitrate;

	f_serial = NULL;
	f_fer = NULL;
	left_limiter_gain = 1.0f;
	right_limiter_gain = 1.0f;
	copyright();
	parsecmdline( argc, argv, &input_filename, &output_filename, &fer_filename, &conf );

	if( conf.FileFormat == F3GP ) {
		// Open3GP(input_filename, 0, &conf); // old
		Open3GP( &tfi, bfi, input_filename, 0, &conf );
	}
	else {
		if( ( f_serial = fopen( input_filename, "rb" ) ) == NULL ) {
			fprintf( stderr, "Error opening the input bitstream file %s.\n", input_filename );
			exit( EXIT_FAILURE );
		}
		ReadHeader( &tfi, bfi, conf.FileFormat, &( conf.extension ), &( conf.mode ), &( conf.st_mode ), &( conf.fscale ), 0, f_serial );
		rewind( f_serial );
	}
	/* In case of bad frame; hard initialization */
	if( ( conf.fs == 16000 || conf.fs == 8000 ) && conf.mode == 14 ) {
		/* Force to wb */
		conf.fscale = 0;
		conf.extension = 0;
	}
	else if( conf.mode == 14 ) {
		/* Force to wb+ if sampling rate desire is different than 16 or 8kHz */
		conf.mode = 2;
		conf.fscale = 96;
		conf.extension = 1;
	}

	/* default sampling rate if undefined */
	if( conf.fs == 0 ) {
		if( conf.extension == 0 )
			conf.fs = 16000;
		else if( conf.fscale == 0 )
			conf.fs = 24000;
		else {
#ifdef FILTER_44kHz
			conf.fs = 44100;
#endif
#ifdef FILTER_48kHz
			conf.fs = 48000;
#endif
		}
	}
	fs_output = conf.fs;

	fac_up = fac_down = 12; /* no oversampling by default */
	frac_down_right = 0;
	frac_down_left = 0;

	if( conf.fscale != 0 ) {
		switch( conf.fs ) {
		case 8000:
			fac_up = 2;
			conf.fs = 48000;
			break;
		case 16000:
			fac_up = 4;
			conf.fs = 48000;
			break;
		case 24000:
			fac_up = 6;
			conf.fs = 48000;
			break;
		case 32000:
			fac_up = 8;
			conf.fs = 48000;
			break;
		case 11025:
			fac_up = 3;
			conf.fs = 44100;
			break;
		case 22050:
			fac_up = 6;
			conf.fs = 44100;
			break;
		}
		set_zero( mem_down_right, 2 * L_FILT_DECIM_FS );
		set_zero( mem_down_left, 2 * L_FILT_DECIM_FS );
	}

	/* read the mode/fscale (1st frame) for initialization */
	select_fs( &conf );

	set_frame_length( &conf, &L_frame );

	num_channels = 2;
	if( ( conf.fscale == 0 ) && ( conf.st_mode < 0 ) || ( conf.extension == 0 ) || ( conf.mono_dec_stereo ) )
		num_channels = 1;

	/* Frame errasure simualtion */
	if( conf.fer_sim ) {
		if( ( f_fer = fopen( fer_filename, "rt" ) ) == NULL ) {
			fprintf( stderr, "Error opening fer file %s.\n", fer_filename );
			exit( 0 );
		}
	}

	/* Open output wav file */
	if( ( f_syn = Wave_fopen( output_filename, "wb", &num_channels, &fs_output, &bitsPerSample, &dataSize ) ) == NULL ) {
		fprintf( stderr, "Error opening output wav file %s.\n", output_filename );
		exit( 0 );
	}

	old_bitrate = get_bitrate( &conf );
	fprintf( stderr, "Decoding @ %6.2fkbps", get_bitrate( &conf ) );

	/* Initialize decoders */
	stAmrwbDec = D_IF_init();
	st_d = malloc( sizeof( Decoder_State_Plus ) );
	init_decoder_amrwb_plus( st_d, (int)num_channels, conf.fscale, 1 );

	fprintf( stderr, "\n --- Running ---\n" );
	/*---------------------------------------------------------------------------*
   * Loop for every analysis/transmission frame.                               *
   *   -New L_FRAME_PLUS data are read. (L_FRAME_PLUS = number of speech data
   *per frame) * -Conversion of the speech data from 16 bit integer to real *
   *   -Call coder_wb to encode the speech.                                    *
   *   -The compressed serial output stream is written to a file.              *
   *   -The synthesis speech is written to a file                              *
   *--------------------------------------------------------------------------*/
	frame = 0;
	tfi = (short)frame % 4;
	old_st_mode = conf.st_mode;

	extension = conf.extension;
	mode = conf.mode;
	st_mode = conf.st_mode;
	fst = conf.fscale;

	while( 1 ) {
		fst_switched = 0;
		fprintf( stderr, "Frames processed: %ld    \r", frame / 4 );

		/* get next frame 4x20ms if AMR-WB+ or 1x20ms if AMR-WB */
		if( conf.FileFormat == F3GP ) {

			if( GetNextFrame3GP( &tfi, bfi, &extension, &mode, &st_mode, &fst, (void *)serial, 0 ) ) {
				break;
			}
			conf.mode = mode;
			conf.st_mode = st_mode;
		}
		else {
			if( !ReadRawFile( &tfi, bfi, &conf, &extension, &mode, &st_mode, &fst, f_serial, (void *)serial ) )
				break;
		}

		/* set the bfi for SIMULATION PURPOSE */
		/* read frame erasures every forth frame for wb+*/
		/* read frame erasures every frame for amr-wb*/
		if( conf.fer_sim ) {
			if( conf.extension > 0 ) {
				for( i = 0; i < 4; i++ ) {
					fscanf( f_fer, "%d", &bfi[i] );
				}
			}
			else {
				fscanf( f_fer, "%d", &bfi[tfi] );
			}
		}

		/* note that mode switching is only allowed within 80ms super-frames */
		if( fst != conf.fscale && !fst_switched ) {
			select_fs( &conf );
			set_frame_length( &conf, &L_frame );
			init_decoder_amrwb_plus( st_d, (int)num_channels, fst, 0 );

			fst_switched = 1;
		}
		/* set new config  */
		if( ( ( extension == 0 ) && ( conf.extension == 1 ) ) || ( ( extension == 1 ) && ( conf.extension == 0 ) ) ) {
			if( ( ( mode >= 0 && mode <= 9 ) || mode == 15 ) && ( conf.extension > 0 ) ) {
				copy_decoder_state( st_d, stAmrwbDec, 1 );
			}
			else if( ( ( mode >= 0 && mode <= 8 ) || mode == 15 ) && ( conf.extension == 0 ) ) {
				copy_decoder_state( st_d, stAmrwbDec, 0 );
			}
			conf.mode = mode;
			conf.extension = extension;
		}

		if( conf.extension > 0 ) {
			nb_samp = decoder_amrwb_plus(
			    (int)conf.mode, (short *)serial, bfi, L_frame, (int)num_channels, channel_right, channel_left, st_d, (int)fst, (int)conf.st_mode, conf.mono_dec_stereo, conf.fscale );

			frame += 4;
		}
		else {
			nb_samp = L_FRAME16k;

			D_IF_decode( stAmrwbDec, (unsigned char *)serial, speech16, bfi[tfi] ? _bad_frame : _good_frame );
			for( j = 0; j < 320; j++ )
				channel_right[j] = speech16[j];
			frame++;
		}

		old_st_mode = conf.st_mode;

		if( num_channels == 2 ) {
			if( conf.limiter_on ) {
				simple_frame_limiter( channel_right, &right_limiter_gain, nb_samp );
				simple_frame_limiter( channel_left, &left_limiter_gain, nb_samp );
			}
			nb_samp_fs = decim_fs( channel_right, nb_samp, channel_right, fac_up, mem_down_right, &frac_down_right );
			nb_samp_fs = decim_fs( channel_left, nb_samp, channel_left, fac_up, mem_down_left, &frac_down_left );

			/* interleave of left and right samples (stereo file format) */
			interleave( channel_right, channel_left, channel_right, nb_samp_fs );
			writ_data( channel_right, 2 * nb_samp_fs, f_syn );
		}
		else {
			if( conf.limiter_on ) {
				simple_frame_limiter( channel_right, &right_limiter_gain, nb_samp );
			}
			nb_samp_fs = decim_fs( channel_right, nb_samp, channel_right, fac_up, mem_down_right, &frac_down_right );
			writ_data( channel_right, nb_samp_fs, f_syn );
		}

		conf.fscale = fst;
		if( fst_switched ) {
			set_frame_length( &conf, &L_frame );
		}

		fflush( f_syn );

		if( fabs( old_bitrate - get_bitrate( &conf ) ) > 0.0001 ) {
			old_bitrate = get_bitrate( &conf );

			if( old_bitrate != 0 )
				fprintf( stderr, "Decoding @ %6.2fkbps\n", old_bitrate );
			else
				/* DTX frame return a bit rate of 0 kbps */
				fprintf( stderr, "Decoding         DTX\n" );
		}
	}

	Wave_fclose( f_syn, bitsPerSample );
	exit( EXIT_SUCCESS );
}
