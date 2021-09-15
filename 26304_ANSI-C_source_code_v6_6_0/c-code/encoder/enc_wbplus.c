#include "../include/amr_plus.h"
#include "../include/wbplus3gplib.h"
#include "../lib_amr/dec_if.h"
#include "../lib_amr/enc_if.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void copyright( void )
{
	fprintf( stderr, "\n" );
	fprintf( stderr, "\n" );
}

static void usage( char *argv )
{
	fprintf( stderr,
	    "Usage: %s -rate <Bit rate> [-mono] | -mi <mode index> [-isf <factor>] [-lc] [-dtx] -ff <3gp/raw> -if <infile.wav> -of <outfile.wb+>\n",
	    argv );
	fprintf( stderr, "\n" );
	fprintf( stderr, "-rate   Bit rate between 6-36 kbps mono or 7-48 kbps stereo \n" );
	fprintf( stderr, "-mono   Force mono encoding \n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "-mi     Mode Index (0..15  -> AMR WB\n                  16..47 -> AMR WB+) (see ts 26.290 Table 25) \n" );
	fprintf( stderr,
	    "-isf    Internal Sampling Frequency (0.5... 1.5, default is 1.0).\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\n" );

	fprintf( stderr, "-lc     low complexity (for AMR-WB+ modes).\n" );
	fprintf( stderr,
	    "-dtx    enables VAD/DTX functionality (for AMR-WB modes).\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "-ff     3gp File Format / raw format\n" );
	fprintf( stderr, "-if     input audio WAV file.\n" );
	fprintf( stderr, "-of     output AMRWB+ 3gp file.\n" );
	fprintf( stderr, "-cf     configuration file\n" );
	fprintf( stderr, "\n" );
}
float get_bitrate( EncoderConfig *conf )
{
	if( conf->fscale != 0 ) {
		return (float)( get_nb_bits( conf->extension, conf->mode, conf->st_mode ) * conf->fscale ) / ( 80.0f * FSCALE_DENOM );
	}
	else {
		return (float)get_nb_bits( conf->extension, conf->mode, conf->st_mode ) / 80.0f;
	}
}

int get_core_mode( float bitrate )
{
	float min_dist = 1e16f;
	int   mode;
	int   i;
	mode = 0;
	for( i = 0; i < 8; i++ ) {
		if( fabs( bitrate - (float)( NBITS_CORE[i] + NBITS_BWE ) / 80.0f ) < min_dist ) {
			min_dist = (float)fabs( bitrate - (float)( NBITS_CORE[i] + NBITS_BWE ) / 80.0f );
			mode = i;
		}
	}
	return mode;
}

int get_stereo_mode( float bitrate )
{
	float min_dist = 1e16f;
	int   mode;
	int   i;

	if( bitrate == 0 )
		return -1;

	mode = 0;
	for( i = 0; i < 16; i++ ) {
		if( fabs( bitrate - (float)( StereoNbits[i] + NBITS_BWE ) / 80.0f ) < min_dist ) {

			min_dist = (float)fabs( bitrate - (float)( StereoNbits[i] + NBITS_BWE ) / 80.0f );
			mode = i;
		}
	}
	return mode;
}

void get_raw_3gp_mode( short *mode, short *st_mode, short raw_3gp_mode, short extension )
{
	short index;
	/* Mono mode only */
	if( raw_3gp_mode <= 8 ) {
		if( extension != 0 ) {
			fprintf( stderr, "-isf is not supported by amr_wb\n" );
			exit( EXIT_FAILURE );
		}

		*mode = raw_3gp_mode;
		*st_mode = -1;
	}
	else if( raw_3gp_mode == 10 ) /* 14m */
	{
		*mode = 2;
		*st_mode = -1;
	}
	else if( raw_3gp_mode == 11 ) /* 18s */
	{
		*mode = 2;
		*st_mode = 6;
	}
	else if( raw_3gp_mode == 12 ) /* 24m*/
	{
		*mode = 7;
		*st_mode = -1;
	}
	else if( raw_3gp_mode == 13 ) /*24s*/
	{
		*mode = 5;
		*st_mode = 7;
	}
	else if( raw_3gp_mode >= 16 && raw_3gp_mode < 24 ) {
		*mode = raw_3gp_mode - 16;
		*st_mode = -1;
	}
	else if( raw_3gp_mode >= 24 && raw_3gp_mode <= 47 ) {
		index = raw_3gp_mode - 24;
		*mode = miMode[2 * index];
		*st_mode = miMode[2 * index + 1];
	}
	else {
		printf( "Invalid Mode Index\n" );
		exit( EXIT_FAILURE );
	}
}

short get_isf_index( short *fscale )
{
	short index, i;
	float dist = 512.0f, ftmp;

	index = 0;

	/* Mono mode only */
	for( i = 0; i < 14; i++ ) {
		ftmp = (float)fabs( *fscale - isfIndex[i] );
		if( ftmp < dist ) {
			dist = ftmp;
			index = i;
		}
	}
	*fscale = isfIndex[index];
	return index;
}

static void parsecmdline( int argc,
    char *                    argv[],
    char **                   input_filename,
    char **                   output_filename,
    char **                   config_filename,
    EncoderConfig *           conf,
    float *                   rate )
{
	int   simple_mode, amr_wb, amr_wbp_carac, mi_mode;
	float srate;
	float mrate;

	if( argc == 1 ) {

		usage( argv[0] );

		exit( EXIT_FAILURE );
	}

	conf->extension = 0;
	conf->allow_dtx = 0;
	conf->use_case_mode = USE_CASE_A;
	conf->fscale = 0;
	conf->mode = -1;
	conf->st_mode = -1;
	conf->FileFormat = F3GP;
	conf->mode_index = -1;
	conf->fscale_index = 0;
	conf->bc = 0;

	simple_mode = 0;
	amr_wb = 0;
	mi_mode = 0;
	amr_wbp_carac = 0;

	mrate = -1;
	srate = -1;
	*rate = -1;
	argc--;
	argv++;

	while( argc > 0 ) {

		if( !strcmp( *argv, "-mi" ) ) {

			if( simple_mode ) {
				fprintf( stderr, "Can't use -rate with -mi\n" );
				exit( EXIT_FAILURE );
			}

			mi_mode = 1;
			argv++;
			argc--;

			conf->mode_index = (short)atoi( *argv );
			if( conf->mode_index < 0 || conf->mode_index > 47 ) {
				fprintf( stderr, "Unknown Mode Index (see TS 26.290 Table 25)\n" );
				exit( EXIT_FAILURE );
			}
			else if( conf->mode_index == 9 || conf->mode_index == 14 || conf->mode_index == 15 ) {
				fprintf( stderr, "Mode Index %d is reserved (see TS 26.290 Table 21)\n", conf->mode_index );
				exit( EXIT_FAILURE );
			}
			else {
				if( ( conf->mode_index >= 0 ) && ( conf->mode_index <= 8 ) ) /* amr_wb modes */
				{
					get_raw_3gp_mode( &( conf->mode ), &( conf->st_mode ), (short)atoi( *argv ), conf->extension = 0 );
					amr_wb = 1;
				}
				else if( ( conf->mode_index >= 10 ) && ( conf->mode_index <= 13 ) ) /* WB+ tested modes */
				{
					get_raw_3gp_mode( &( conf->mode ), &( conf->st_mode ), (short)atoi( *argv ), conf->extension = 0 );
					conf->extension = 1;
					amr_wbp_carac = 1;
				}
				else {
					conf->extension = 1;
					get_raw_3gp_mode( &( conf->mode ), &( conf->st_mode ), conf->mode_index, conf->extension );

					if( conf->fscale == 0 ) {
						conf->fscale = FSCALE_DENOM;
						conf->fscale_index = 8;
					}
				}
			}
		}
		else if( !strcmp( *argv, "-isf" ) ) {
			if( simple_mode ) {
				fprintf( stderr, "Can't use -rate with -isf\n" );
				exit( EXIT_FAILURE );
			}
			if( amr_wb ) {
				fprintf( stderr, "-isf is not supported by amr_wb\n" );
				exit( EXIT_FAILURE );
			}
			argv++;
			argc--;

			mi_mode = 1; /* -isf is only allow with -mi */
			conf->extension = 1;
			if( ( atof( *argv ) >= 0.5 ) && ( atof( *argv ) <= 1.5 ) ) {
				conf->fscale = (short)( ( atof( *argv ) * FSCALE_DENOM ) + 0.5f );
				/* force scale to be an even number */
				conf->fscale = ( conf->fscale >> 1 ) << 1;
				/* limit the scale factor */
				if( conf->fscale > FAC_FSCALE_MAX ) {
					conf->fscale = FAC_FSCALE_MAX;
				}
				if( conf->fscale < FAC_FSCALE_MIN ) {
					conf->fscale = FAC_FSCALE_MIN;
				}
				conf->fscale_index = get_isf_index( &( conf->fscale ) ); /* Use "fscale from index" */
			}
			else {
				fprintf( stderr, "Unknown Inernal Sampling Frequency factor\n" );
				exit( EXIT_FAILURE );
			}
		}
		else if( !strcmp( *argv, "-rate" ) ) {
			if( mi_mode ) {
				fprintf( stderr, "Can't use -rate with -mi or -isf \n" );
				exit( EXIT_FAILURE );
			}
			argv++;
			argc--;
			simple_mode = 1;
			*rate = (float)atof( *argv );
			conf->extension = 1;
			if( *rate < 6.0 || *rate > 48.0 ) {
				fprintf( stderr, "Minimum rate is 6.0kbps and maximum rate is 48.0 kbps\n" );
				exit( EXIT_FAILURE );
			}
		}
		else if( !strcmp( *argv, "-mono" ) ) {
			conf->st_mode = -2; /* indicate mono is forced */
			if( *rate > 36.0 ) {
				fprintf( stderr, "Maximum mono rate is 36.0 kbps\n" );
				exit( EXIT_FAILURE );
			}
		}
		else if( !strcmp( *argv, "-lc" ) ) {
			conf->use_case_mode = USE_CASE_B;
		}
		else if( !strcmp( *argv, "-dtx" ) ) {
			conf->allow_dtx = 1;
		}
		else if( !strcmp( *argv, "-bc" ) ) {
			conf->bc = 1;
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
		else if( !strcmp( *argv, "-cf" ) ) {
			argv++;
			argc--;
			*config_filename = *argv;
		}
		else if( !strcmp( *argv, "-ff" ) ) {
			argv++;
			argc--;
			if( !strcmp( *argv, "raw" ) ) {
				conf->FileFormat = FRAW;
			}
			else {
				conf->FileFormat = F3GP;
			}
		}
		else {
			fprintf( stderr, "Unknown option %s\n", *argv );
			exit( EXIT_FAILURE );
		}
		argv++;
		argc--;
	}
	if( amr_wbp_carac && conf->fscale != 0 ) {
		fprintf( stderr, "-isf is not supported with AMR WB caracterized modes\n" );
		exit( EXIT_FAILURE );
	}
	if( conf->st_mode == -2 && simple_mode != 1 ) {
		fprintf( stderr, "Choose right Mode Index to encode mono File\n-mono option is only supported with -rate\n" );
		exit( EXIT_FAILURE );
	}
}

static void set_frame_length( int samplingRate,
    int                           fscale,
    int *                         L_frame,
    int *                         L_next,
    int *                         L_next_st )
{
	if( fscale != 0 ) {
		switch( samplingRate ) {
#ifdef FILTER_48kHz
		case 8000:
		case 16000:
		case 24000:
		case 32000:
		case 48000:
			*L_frame = 2 * L_FRAME48k;
			break;
#endif
#ifdef FILTER_44kHz
		case 11025:
		case 22050:
		case 44100:
			*L_frame = 2 * L_FRAME44k;
			break;
#endif
		default:
			fprintf( stderr, "error in sampling frequency.  choice of filter are: \n" );
#ifdef FILTER_44kHz
			fprintf( stderr, "     11, 22, 44 kHz \n" );
#endif
#ifdef FILTER_48kHz
			fprintf( stderr, "     8, 16, 24, 32, 48 kHz \n" );
#endif
			exit( EXIT_FAILURE );
			break;
		}
	}
	else {
		switch( samplingRate ) {
		case 8000:
			*L_frame = L_FRAME8k;
			*L_next = L_NEXT8k;
			*L_next_st = L_NEXT_ST8k;
			break;
		case 16000:
			*L_frame = L_FRAME16kPLUS;
			*L_next = L_NEXT16k;
			*L_next_st = L_NEXT_ST16k;
			break;
		case 24000:
			*L_frame = L_FRAME24k;
			*L_next = L_NEXT24k;
			*L_next_st = L_NEXT_ST24k;
			break;
		default:
			fprintf( stderr, "error in sampling freq: without fsratio(isf) only 8, 16 or 24 kHz are allowed\n" );
			exit( EXIT_FAILURE );
			break;
		}
	}
}

static void GetRate( EncoderConfig *conf, float rate, const short *TableRate, short lenght )
{
	short index, i;
	float dist = 512.0f, ftmp;

	index = 0;
	/* Mono mode only */
	for( i = 0; i < lenght; i += 3 ) {
		ftmp = (float)fabs( rate * 2 - TableRate[i] );
		if( ftmp < dist ) {
			dist = ftmp;
			index = i;
		}
	}
	conf->mode_index = TableRate[index + 1];
	conf->fscale_index = TableRate[index + 2];
	get_raw_3gp_mode( &conf->mode, &conf->st_mode, conf->mode_index, conf->extension );
	conf->fscale = isfIndex[conf->fscale_index];
}

static void moveAndRound( float *in, short *out, int n )
{
	int   i;
	float temp;

	for( i = 0; i < n; i++ ) {
		temp = *in++;
		if( temp >= 0.0 )
			temp += 0.5;
		else
			temp -= 0.5;
		if( temp > 32767.0 )
			temp = 32767.0;
		if( temp < -32767.0 )
			temp = -32767.0;
		*out++ = (short)temp;
	}
}
static void deinterleave( float *buf, float *left, float *right, int length )
{
	int i;

	for( i = 0; i < length; i++ ) {
		left[i] = buf[i * 2];
		right[i] = buf[( i * 2 ) + 1];
	}
}

int get_config( FILE *fp, float t[] )
{
	int OK = 0;
	if( !fp || feof( fp ) )
		return 0;

	/*
    read from config file the following items
    time  extension  mode_index fscale
    time (in seconds) must always be above 0.0, 0.000001 is OK
  */
	while( !OK && !feof( fp ) ) {
		char s[100], *sp;
		int  ix = 0;
		t[0] = 0.0;
		fgets( s, 99, fp );
		sp = strtok( s, " \t" );
		while( sp && ix < 4 ) {
			t[ix++] = (float)atof( sp );

			sp = strtok( 0, " \t" );
		}
		if( t[0] != 0.0 ) {
			return 1;
		}
	}
	return 0;
}
void close_wbp( Coder_State_Plus *st, Word16 UseCaseB )
{

	if( st->stClass != NULL && UseCaseB > 0 ) {
		free( st->stClass );
		st->stClass = NULL;
	}
	if( st->vadSt != NULL && UseCaseB > 0 )
		wb_vad_exit( &st->vadSt );

	if( st != NULL ) {
		free( st );
		st = NULL;
	}
}
void main( int argc, char *argv[] )
{
	FILE *                f_speech; /* File of speech data                   */
	FILE *                f_serial; /* File of serial bits for transmission  */
	short                 speech16[L_FRAME_FSMAX * 2];
	short                 serial[NBITS_MAX];      /* serial parameters.                    */
	unsigned char         serialAmrwb[NBITS_MAX]; /* serial parameters.                    */
	float                 channel_right[4 * L_FRAME_FSMAX];
	float                 channel_left[2 * L_FRAME_FSMAX];
	float                 mem_up_right[2 * L_FILT_OVER_FS], mem_up_left[2 * L_FILT_OVER_FS];
	int                   frac_up_right, frac_up_left;
	int                   fac_up, fac_down, nb_samp_fs;
	Coder_State_Plus *    st = NULL;
	short                 numOfChannels, bitsPerSample;
	void *                stAmrwbEnc = NULL;
	int                   nb_bits = 0;
	int                   i, lg, L_frame;
	long                  frame, samplingRate, dataSize;
	int                   L_next, L_next_st;
	int                   nb_samp, nb_hold;
	short /*codec_mode,*/ mode, extension;
	short                 serial_size = 0;
	char *                input_filename;
	char *                output_filename;
	int                   st_mode, old_st_mode;
	EncoderConfig         conf;
	short                 fst;
	float                 rate;
	char                  FileFormatType[25];

	FILE *f_config = 0;
	float rec_time = 0.0;
	float config_file_time = 0.0;
	float t[4]; /* time,conf.extension,conf.mode_index,conf.fscale */
	float old_bitrate;
	char *config_filename;

	int Processed_sample = 0;
	int nb_sample_to_process = 0;

	f_serial = NULL;
	/* Initializations */
	input_filename = NULL;
	output_filename = NULL;
	config_filename = NULL;

	copyright();
	parsecmdline( argc, argv, &input_filename, &output_filename, &config_filename, &conf, &rate );

	/* Open input wave file */
	if( ( f_speech = Wave_fopen( input_filename, "rb", &numOfChannels, &samplingRate, &bitsPerSample, &dataSize ) ) == NULL ) {
		fprintf( stderr, "Error opening the input file %s.\n", input_filename );
		exit( EXIT_FAILURE );
	}
	/* Simple interface */
	if( rate != -1.0 ) {
		if( conf.st_mode == -2 || numOfChannels == 1 ) {
			GetRate( &conf, rate, MonoRate, 3 * 18 );
		}
		else {
			GetRate( &conf, rate, StereoRate, 3 * 26 );
		}
	}
	/* test if it is stereo input */
	if( ( conf.st_mode >= 0 ) && ( numOfChannels != 2 ) ) {
		fprintf( stderr, "Input file %s must be stereo\n", input_filename );
		exit( EXIT_FAILURE );
	}
	/* test if it is 16 bits PCM */
	if( bitsPerSample != 16 ) {
		fprintf( stderr, "Input file %s must be 16 bits encoded\n", input_filename );
		exit( EXIT_FAILURE );
	}

	if( conf.FileFormat == F3GP ) {
		strcpy( FileFormatType, "3gp File Format" );
		if( conf.bc == 1 ) {
			/* create backward compatible file */
			Create3GPAMRWB();
		}
		else {
			Create3GPAMRWBPlus();
		}
	}
	else /* raw data */
	{
		strcpy( FileFormatType, "Raw File Format" );

		/* Open the output bitstream file */
		if( ( f_serial = fopen( output_filename, "wb" ) ) == NULL ) {
			fprintf( stderr, "Error opening output bitstream file %s.\n", output_filename );
			exit( EXIT_FAILURE );
		}
	}

	if( ( conf.extension == 0 ) && ( samplingRate != 16000 ) ) {
		fprintf( stderr, "AMR-WB work only at 16kHz\n" );
		exit( EXIT_FAILURE );
	}

	fac_up = fac_down = 12; /* no oversampling by default */
	frac_up_right = 0;
	frac_up_left = 0;

	if( conf.fscale != 0 ) {
		switch( samplingRate ) {
		case 8000:
			fac_down = 2;
			samplingRate = 48000;
			break;
		case 16000:
			fac_down = 4;
			samplingRate = 48000;
			break;
		case 24000:
			fac_down = 6;
			samplingRate = 48000;
			break;
		case 32000:
			fac_down = 8;
			samplingRate = 48000;
			break;
		case 11025:
			fac_down = 3;
			samplingRate = 44100;
			break;
		case 22050:
			fac_down = 6;
			samplingRate = 44100;
			break;
		}
		set_zero( mem_up_right, 2 * L_FILT_OVER_FS );
		set_zero( mem_up_left, 2 * L_FILT_OVER_FS );
	}

	/* Set default buffer lengths */
	set_frame_length( samplingRate, conf.fscale, &L_frame, &L_next, &L_next_st );
	st = malloc( sizeof( Coder_State_Plus ) );
	memset( serial, 0x42, NBITS_MAX * sizeof( short ) );

	old_bitrate = get_bitrate( &conf );
	fprintf( stderr, "%s\nEncoding @ %6.2fkbps", FileFormatType, get_bitrate( &conf ) );

	init_coder_amrwb_plus( st, (int)numOfChannels, conf.fscale, conf.use_case_mode, 1 );
	stAmrwbEnc = E_IF_init();

	/* number of sample per channel to read from file */
	nb_samp_fs = ( ( L_frame * fac_down ) + frac_up_right ) / fac_up;

	lg = read_data( f_speech, channel_right, ( numOfChannels * nb_samp_fs ) );
	if( lg != ( numOfChannels * nb_samp_fs ) ) {
		printf( "Error: file too short!\n" );
		exit( EXIT_FAILURE );
	}
	if( numOfChannels == 2 ) {
		deinterleave( channel_right, channel_left, channel_right, nb_samp_fs );

		over_fs( channel_left, channel_left, L_frame, fac_down, mem_up_left, &frac_up_left );
	}
	over_fs( channel_right, channel_right, L_frame, fac_down, mem_up_right, &frac_up_right );

	if( conf.extension > 0 ) {
		nb_samp = coder_amrwb_plus_first( channel_right, channel_left, numOfChannels, L_frame, ( numOfChannels == 1 ) ? L_next : L_next_st, conf.fscale, st );
	}
	else {
		moveAndRound( channel_right, speech16, 320 );
		E_IF_encode_first( stAmrwbEnc, speech16 );
		nb_samp = 320;
	}
	fprintf( stderr, "\n --- Running ---\n" );

	/*---------------------------------------------------------------------------*
	* Loop for every analysis/transmission frame.                               *
	*   -New L_FRAME_PLUS data are read. (L_FRAME_PLUS = number of speech data per frame) *
	*   -Conversion of the speech data from 16 bit integer to real              *
	*   -Call coder_wb to encode the speech.                                    *
	*   -The compressed serial output stream is written to a file.              *
	*   -The synthesis speech is written to a file                              *
	*--------------------------------------------------------------------------*/
	mode = conf.mode;
	extension = conf.extension;
	old_st_mode = st_mode = conf.st_mode;

	frame = 0;
	fst = conf.fscale;
	if( config_filename != 0 ) {
		if( fac_down != 12 ) {
			fprintf( stderr,
			    "Must have 16kHz (to switch amr-wb -> wb+ (mode 10@13))\n"
			    "or 48kHz (wb+) input\n\n" );
			exit( EXIT_FAILURE );
		}
		if( conf.mode_index < 16 && samplingRate != 16000 ) {
			fprintf( stderr, "Must have a 16kHz input file (to switch amr-wb -> wb+ (mode 10@13))\n" );
			exit( EXIT_FAILURE );
		}
		if( numOfChannels == 1 ) {
			fprintf( stderr, "** Warning **\nDo not switch to stereo... \nyou have mono input\n\n" );
		}
		f_config = fopen( config_filename, "r" );
		if( !f_config ) {
			fprintf( stderr, "Error opening config file %s.\n", config_filename );
			exit( EXIT_FAILURE );
		}
		while( !get_config( f_config, t ) && !feof( f_config ) ) {
			printf( "%2.3f \n", t[0] );
		}
		config_file_time = t[0];
		rec_time = 0;
	}
	nb_sample_to_process = dataSize * numOfChannels - lg;
	Processed_sample += nb_samp;

	while( Processed_sample < nb_sample_to_process || frame == 0 ) {
		fprintf( stderr, " Frames processed: %ld    \r", frame );
		frame++;

		if( f_config ) {
			if( rec_time >= config_file_time && ( t[0] != -1.0 ) && !feof( f_config ) ) {
				short tmp_mode_index;
				extension = (short)( t[1] + 0.01 );
				tmp_mode_index = (short)( t[2] + 0.01 );
				if( extension == 0 ) {
					mode = tmp_mode_index;
					st_mode = -1;
					conf.mode = tmp_mode_index;
					conf.st_mode = -1;
					conf.mode_index = tmp_mode_index;
				}
				else {
					conf.mode_index = tmp_mode_index;
					get_raw_3gp_mode( &( conf.mode ), &( conf.st_mode ), conf.mode_index, extension );
				}
				fst = (short)( t[3] * FSCALE_DENOM + 0.5 );

				get_isf_index( &fst ); /* Use "fscale from index" */

				while( !get_config( f_config, t ) && !feof( f_config ) ) {
				}
				config_file_time = t[0];
			}
		}
		if( fst != conf.fscale ) {
			conf.fscale = fst;
			init_coder_amrwb_plus( st, (int)numOfChannels, conf.fscale, conf.use_case_mode, 0 );
			set_frame_length( samplingRate, conf.fscale, &L_frame, &L_next, &L_next_st );
			conf.fscale_index = get_isf_index( &( conf.fscale ) ); /* Use "fscale from index" */
		}

		nb_hold = L_frame - nb_samp;

		mvr2r( channel_right + nb_samp, channel_right, nb_hold );
		mvr2r( channel_left + nb_samp, channel_left, nb_hold );

		/* number of sample per channel to read from file */
		nb_samp_fs = ( ( nb_samp * fac_down ) + frac_up_right ) / fac_up;

		lg = read_data( f_speech, channel_right + nb_hold, ( numOfChannels * nb_samp_fs ) );

		if( numOfChannels == 2 ) {
			deinterleave( channel_right + nb_hold, channel_left + nb_hold, channel_right + nb_hold, nb_samp_fs );

			over_fs( channel_left + nb_hold, channel_left + nb_hold, nb_samp, fac_down, mem_up_left, &frac_up_left );
		}
		over_fs( channel_right + nb_hold, channel_right + nb_hold, nb_samp, fac_down, mem_up_right, &frac_up_right );

		if( ( ( extension == 0 ) && ( conf.extension == 1 ) )
		    || ( ( extension == 1 ) && ( conf.extension == 0 ) ) ) {
			if( ( ( mode >= 0 && mode <= 9 ) || mode == 15 ) && ( conf.extension > 0 ) ) {
				copy_coder_state( st, stAmrwbEnc, 1, conf.use_case_mode );
			}
			else if( ( mode >= 0 && mode <= 8 ) && ( conf.extension == 0 ) ) {
				copy_coder_state( st, stAmrwbEnc, 0, conf.use_case_mode );
			}
			/* conf.mode = mode; */
			conf.extension = extension;
		}
		if( conf.extension > 0 ) {
			/* update needed if mode changes */
			nb_bits = get_nb_bits( conf.extension, conf.mode, conf.st_mode );
			if( numOfChannels == 2 ) {
				nb_samp = coder_amrwb_plus_stereo( channel_right, channel_left, conf.mode, L_frame, serial, st, conf.use_case_mode, conf.fscale, conf.st_mode );
			}
			else {
				nb_samp = coder_amrwb_plus_mono( channel_right,
				    conf.mode,
				    L_frame,
				    serial,
				    st,
				    conf.use_case_mode,
				    conf.fscale );
			}
			old_st_mode = conf.st_mode;

			if( conf.FileFormat == F3GP ) {
				WriteSamplesAMRWBPlus( conf, serial, nb_bits );
			}
			else {
				for( i = 0; i < 4; i++ ) {
					WriteHeader( conf, (short)nb_bits, (short)i, f_serial );
					WriteBitstreamPlus( conf, (short)nb_bits, (short)i, serial, f_serial );
				}
			}
		}
		else {
			for( i = 0; i < 4; i++ ) {
				moveAndRound( &channel_right[i * 320], speech16, 320 );
				serial_size = (short)E_IF_encode( stAmrwbEnc, (Word16)conf.mode, speech16, serialAmrwb, conf.allow_dtx );
				if( conf.FileFormat == F3GP ) {
					WriteSamplesAMRWBPlus( conf, serialAmrwb, serial_size );
				}
				else {
					WriteHeader( conf, (short)serial_size, (short)i, f_serial );
					WriteBitstream( conf, (short)serial_size, (short)i, serialAmrwb, f_serial );
				}
			}
			nb_samp = 4 * L_FRAME16k;
		}
		if( fabs( old_bitrate - get_bitrate( &conf ) ) > 0.00001 ) {
			old_bitrate = get_bitrate( &conf );
			fprintf( stderr, "Rectime: %2.3f Encoding @ %6.2fkbps\n", rec_time, old_bitrate );
		}
		rec_time += nb_samp / ( (float)( samplingRate ) );
		Processed_sample += ( nb_samp_fs * numOfChannels );
	}

	if( conf.FileFormat == F3GP ) {
		Close3GP( output_filename );
	}
	else {
		fclose( f_serial );
	}
	Wave_fclose( f_speech, bitsPerSample );

	if( stAmrwbEnc != NULL )
		E_IF_exit( stAmrwbEnc );

	close_wbp( st, conf.use_case_mode );

	exit( EXIT_SUCCESS );
}
