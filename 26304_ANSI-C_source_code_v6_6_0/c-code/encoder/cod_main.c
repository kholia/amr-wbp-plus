#include "../include/amr_plus.h"
#include "../include/s_util.h"
#include "../lib_amr/dec_if.h"
#include "../lib_amr/enc_if.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
 * Funtion  init_coder_amr_plus                                    *
 *          ~~~~~~~~~~~~~~~~~~~                                    *
 *   - Allocate memory for static variables.                       *
 *   - Initialization of variables for the coder section.          *
 *   - compute lag window and LP analysis window (if not done)     *
 *-----------------------------------------------------------------*/
void init_coder_amrwb_plus( Coder_State_Plus *st, int num_chan, int fscale, short use_case_mode, short full_reset )
{
	if( full_reset && use_case_mode == 1 ) {
		/* initialize memories (ClassB part) */
		st->vad_hist = 0;
		wb_vad_init( &( st->vadSt ) );
		st->stClass = malloc( sizeof( NCLASSDATA ) );
		initClassifyExcitation( st->stClass );
	}
	/* initialize bwe codebooks */
	if( fscale == 0 ) {
		st->left.mean_isf_hf = (float *)mean_isf_hf_low_rate;
		st->left.dico1_isf_hf = (float *)dico1_isf_hf_low_rate;
		st->right.mean_isf_hf = (float *)mean_isf_hf_low_rate;
		st->right.dico1_isf_hf = (float *)dico1_isf_hf_low_rate;
	}
	else {
		st->left.mean_isf_hf = (float *)mean_isf_hf_12k8;
		st->left.dico1_isf_hf = (float *)dico1_isf_hf_12k8;
		st->right.mean_isf_hf = (float *)mean_isf_hf_12k8;
		st->right.dico1_isf_hf = (float *)dico1_isf_hf_12k8;
	}
	/* initialize memories (stereo part) */
	st->right.decim_frac = 0;
	st->left.decim_frac = 0;
	if( full_reset ) {
		/* initialize memories (stereo part) */
		set_zero( st->right.mem_decim_hf, 2 * L_FILT24k );
		set_zero( st->right.old_speech_hf, L_OLD_SPEECH_ST );
		set_zero( st->right.mem_decim, L_MEM_DECIM_SPLIT );
		st->right.mem_preemph = 0.0;
		init_coder_hf( &( st->right ) );
		set_zero( st->left.mem_decim, L_MEM_DECIM_SPLIT );
		set_zero( st->left.mem_decim_hf, 2 * L_FILT24k );
		set_zero( st->left.old_speech_hf, L_OLD_SPEECH_ST );
		st->left.mem_preemph = 0.0;
		init_coder_hf( &( st->left ) );
		/* initialize memories (mono part) */
		set_zero( st->old_speech, L_OLD_SPEECH_ST );
		set_zero( st->old_synth, M );
		set_zero( st->left.mem_sig_in, 4 );
		set_zero( st->right.mem_sig_in, 4 );
		/* init band splitting memories */
		set_zero( st->old_chan, L_OLD_SPEECH_ST );
		set_zero( st->old_chan_2k, L_OLD_SPEECH_2k );
		set_zero( st->old_speech_2k, L_OLD_SPEECH_2k );
		set_zero( st->old_chan_hi, L_OLD_SPEECH_hi );
		set_zero( st->old_speech_hi, L_OLD_SPEECH_hi );
		init_coder_lf( st );

		if( num_chan == 2 ) {
			init_coder_stereo_x( st );
		}
		st->SwitchFlagPlusToWB = 0;
		set_zero( st->mem_gain_code, 4 );
		st->prev_mod = 0;
	}
	else {
	}
	/* Initialize the LP analysis window */
	if( fscale <= FSCALE_DENOM ) {
		cos_window( st->window, L_WINDOW / 2, L_WINDOW / 2 );
	}
	else {
		cos_window( st->window, L_WINDOW_HIGH_RATE / 2, L_WINDOW_HIGH_RATE / 2 );
	}
	return;
}
/*-----------------------------------------------------------------*
 * Funtion  coder_amrwb_plus_first                                 *
 *          ~~~~~~~~~~~~~~~~                                       *
 *   - Fill lookahead buffers                                      *
 *                                                                 *
 *-----------------------------------------------------------------*/
int coder_amrwb_plus_first(            /* output: number of sample processed */
    float             channel_right[], /* input: used on mono and stereo       */
    float             channel_left[],  /* input: used on stereo only           */
    int               n_channel,       /* input: 1 or 2 (mono/stereo)              */
    int               L_frame,         /* input: frame size                        */
    int               L_next,          /* input: lookahead                         */
    int               fscale,
    Coder_State_Plus *st /* i/o : coder memory state                 */
)
{
	float  buffer[L_FRAME_FSMAX];
	float  old_speech[L_TOTAL_ST];
	float *new_speech = old_speech + L_OLD_SPEECH_ST;
	float  sig_right[L_FRAME_PLUS];
	float  old_mono_2k[L_TOTAL_ST_2k];
	float *new_mono_2k = old_mono_2k + L_OLD_SPEECH_2k;
	float  old_mono_hi[L_TOTAL_ST_hi];
	float *new_mono_hi = old_mono_hi + L_OLD_SPEECH_hi;
	int    nb_samp, fac_fs;

	int WorkLen, fac_up, fac_down;

	if( ( L_frame - L_FRAME32k ) == 0 ) {
		fac_fs = FSCALE_DENOM * 3 / 2;
	}
	else {
		fac_fs = fscale;
	}
	/* 48k setting*/
	fac_up = ( fac_fs << 3 );
	fac_down = 180 * 8;

	if( fscale != 0 ) {

#ifdef FILTER_44kHz
		if( ( L_frame - ( 2 * L_FRAME44k ) ) == 0 ) {
			fac_up = ( fac_fs << 3 );
			fac_down = 3 * 441;
		}
#endif
		if( n_channel == 2 ) {
			/*L_frame = ((L_frame_int*fac_down)+(*frac_mem))/fac_up;*/
			WorkLen = ( ( L_NEXT_ST * 2 * fac_down ) + ( st->right.decim_frac ) ) / fac_up;
		}
		else {
			/*L_frame = ((L_frame_int*fac_down)+(*frac_mem))/fac_up;*/
			WorkLen = ( ( L_NEXT * 2 * fac_down ) + ( st->right.decim_frac ) ) / fac_up;
		}
	}
	else {
		WorkLen = L_next;
	}

	/*-----------------------------------------------------------------*
  * MONO/STEREO signal downsampling (codec working at 6.4kHz)       *
  * - decimate signal to fs=12.8kHz                                 *
  * - Perform 50Hz HP filtering of signal at fs=12.8kHz.            *
  * - Perform fixed preemphasis through 1 - g z^-1                  *
  * - Mix left and right channels into sum and difference signals   *
  *-----------------------------------------------------------------*/
	set_zero( buffer, L_FRAME_FSMAX );
	set_zero( old_speech, L_OLD_SPEECH_ST );
	if( n_channel == 2 ) {
		set_zero( old_mono_2k, L_OLD_SPEECH_2k );
		set_zero( old_mono_hi, L_OLD_SPEECH_hi );
		if( fscale == 0 ) {
			/* copy memory into working space */
			mvr2r( channel_left, buffer + L_FRAME_FSMAX - L_next, L_next );
			decim_12k8( buffer + L_FRAME_FSMAX - L_frame, L_frame, new_speech, st->left.mem_decim, 0 );
			/* copy memory into working space */
			mvr2r( channel_right, buffer + L_FRAME_FSMAX - L_next, L_next );
			decim_12k8( buffer + L_FRAME_FSMAX - L_frame, L_frame, sig_right, st->right.mem_decim, 0 );
			nb_samp = L_next;
		}
		else {
			set_zero( new_speech, L_FRAME_PLUS - L_NEXT_ST );
			set_zero( sig_right, L_FRAME_PLUS - L_NEXT_ST );

			decim_split_12k8( channel_left, L_frame, new_speech + L_FRAME_PLUS - L_NEXT_ST, channel_left, L_NEXT_ST, fac_fs, fac_up, fac_down, WorkLen, st->left.mem_decim, &( st->left.decim_frac ) );

			nb_samp = decim_split_12k8( channel_right, L_frame, sig_right + L_FRAME_PLUS - L_NEXT_ST, channel_right, L_NEXT_ST, fac_fs, fac_up, fac_down, WorkLen, st->right.mem_decim, &( st->right.decim_frac ) );
		}

		hp50_12k8( new_speech, L_FRAME_PLUS, st->left.mem_sig_in, fscale );
		hp50_12k8( sig_right, L_FRAME_PLUS, st->right.mem_sig_in, fscale );
		/* parametric stereo : mix left and right channels */
		mix_ch( new_speech, sig_right, new_speech, L_FRAME_PLUS, 1.0f, 1.0f );
		/* do the lo,hi band-splitting on the mono signal */
		band_split_taligned_2k( new_speech, new_mono_2k, new_mono_hi, L_FRAME_PLUS );
		/* copy working space into memory */
		mvr2r( old_speech + L_FRAME_PLUS, st->old_speech, L_OLD_SPEECH_ST );
		mvr2r( old_mono_2k + L_FRAME_2k, st->old_speech_2k, L_OLD_SPEECH_2k );
		mvr2r( old_mono_hi + L_FRAME_PLUS, st->old_speech_hi, L_OLD_SPEECH_hi );
	}
	else {
		if( fscale == 0 ) {
			/* copy memory into working space */
			mvr2r( channel_right, buffer + L_FRAME_FSMAX - L_next, L_next );
			decim_12k8( buffer + L_FRAME_FSMAX - L_frame, L_frame, new_speech, st->right.mem_decim, 0 );
			nb_samp = L_next;
		}
		else {
			set_zero( new_speech, L_FRAME_PLUS - L_NEXT );
			/* decimation and band split (HF temporary into channel_right) */
			nb_samp = decim_split_12k8( channel_right, L_frame, new_speech + L_FRAME_PLUS - L_NEXT, channel_right, L_NEXT, fac_fs, fac_up, fac_down, WorkLen, st->right.mem_decim, &( st->right.decim_frac ) );
		}
		hp50_12k8( new_speech, L_FRAME_PLUS, st->right.mem_sig_in, fscale );
		/* copy working space into memory */
		mvr2r( old_speech + L_FRAME_PLUS, st->old_speech, L_OLD_SPEECH_ST );
	}
	/* Apply preemphasis (for core codec only */
	E_UTIL_f_preemph( new_speech, PREEMPH_FAC, L_FRAME_PLUS, &( st->right.mem_preemph ) );
	/* update lower band memory for next frame */
	mvr2r( old_speech + L_FRAME_PLUS, st->old_speech_pe, L_OLD_SPEECH_ST );
	if( n_channel == 2 ) {
		mvr2r( sig_right, new_speech, L_FRAME_PLUS );
		/* do the lo,hi band-splitting on the mono signal */
		band_split_taligned_2k( new_speech, new_mono_2k, new_mono_hi, L_FRAME_PLUS );
		/* copy working space into memory */
		mvr2r( old_speech + L_FRAME_PLUS, st->old_chan, L_OLD_SPEECH_ST );
		mvr2r( old_mono_2k + L_FRAME_2k, st->old_chan_2k, L_OLD_SPEECH_2k );
		mvr2r( old_mono_hi + L_FRAME_PLUS, st->old_chan_hi, L_OLD_SPEECH_hi );
	}
	if( L_frame > L_FRAME8k ) {
		/* prepare buffers for MONO/STEREO Bandwidth extension */
		if( fscale == 0 ) {
			/* copy memory into working space */
			mvr2r( channel_right, buffer + L_FRAME_FSMAX - L_next, L_next );
			decim_12k8( buffer + L_FRAME_FSMAX - L_frame, L_frame, new_speech, st->right.mem_decim_hf, ( fscale == 0 ) ? 1 : 2 );
		}
		else {
			if( n_channel == 2 ) {
				/* right HF was stored into channel_right */
				mvr2r( channel_right, new_speech + L_FRAME_PLUS - L_NEXT_ST, L_NEXT_ST );
			}
			else {
				/* right HF was stored into channel_right */
				mvr2r( channel_right, new_speech + L_FRAME_PLUS - L_NEXT, L_NEXT );
			}
		}
		mvr2r( old_speech + L_FRAME_PLUS, st->right.old_speech_hf, L_OLD_SPEECH_ST );
		if( n_channel == 2 ) {
			if( fscale == 0 ) {
				/* copy memory into working space */
				mvr2r( channel_left, buffer + L_FRAME_FSMAX - L_next, L_next );
				decim_12k8( buffer + L_FRAME_FSMAX - L_frame, L_frame, new_speech, st->left.mem_decim_hf, ( fscale == 0 ) ? 1 : 2 );
			}
			else {
				/* left HF was stored into channel_left */
				mvr2r( channel_left, new_speech + L_FRAME_PLUS - L_NEXT_ST, L_NEXT_ST );
			}
			mvr2r( old_speech + L_FRAME_PLUS, st->left.old_speech_hf, L_OLD_SPEECH_ST );
		}
	}
	return ( nb_samp );
}
/*-----------------------------------------------------------------*
 * Funtion  coder_amrwb_plus_stereo                                *
 *          ~~~~~~~~~~~~~~~~                                       *
 *   - Main stereo coder routine.                                  *
 *                                                                 *
 *-----------------------------------------------------------------*/
int coder_amrwb_plus_stereo(           /* output: number of sample processed */
    float             channel_right[], /* input: used on mono and stereo       */
    float             channel_left[],  /* input: used on stereo only           */
    int               codec_mode,      /* input: AMR-WB+ mode (see cnst.h)         */
    int               L_frame,         /* input: 80ms frame size                   */
    short             serial[],        /* output: serial parameters                */
    Coder_State_Plus *st,              /* i/o : coder memory state                 */
    short             use_case_mode,
    int               fscale,
    int               StbrMode )
{
	/* LPC coefficients of lower frequency */
	float AqLF[( NB_SUBFR + 1 ) * ( M + 1 )];
	int   param[NB_DIV * NPRM_DIV];
	int   prm_stereo[MAX_NPRM_STEREO_DIV * NB_DIV]; /* see cnst.h */
	int   prm_hf_left[NB_DIV * NPRM_BWE_DIV];
	int   prm_hf_right[NB_DIV * NPRM_BWE_DIV];
	/* vector working at fs=12.8kHz */
	float  sig_left[L_FRAME_PLUS];
	float  sig_right[L_FRAME_PLUS];
	float  old_speech[L_TOTAL_ST];
	float *speech, *new_speech;
	float  old_synth[M + L_FRAME_PLUS];
	float *synth;
	float *syn_hf = old_synth;
	/* Scalars */
	int i, k, nbits_pack;
	int mod[NB_DIV];
	/* LTP parameters for high band */
	float ol_gain[NB_DIV];
	int   T_out[NB_DIV];
	float p_out[NB_DIV];
	int   nb_samp, fac_fs;
	/*ClassB parameters*/
	short excType[4];
	int   WorkLen, fac_up, fac_down;

	/* decimation and band split (HF temporary into channel_right) */
	if( ( L_frame - L_FRAME32k ) == 0 ) {
		fac_fs = FSCALE_DENOM * 3 / 2;
	}
	else {
		fac_fs = fscale;
	}

	/* 48k setting*/
	fac_up = ( fac_fs << 3 );
	fac_down = 180 * 8;

	if( fscale != 0 ) {

#ifdef FILTER_44kHz
		if( ( L_frame - ( 2 * L_FRAME44k ) ) == 0 ) {
			fac_up = ( fac_fs << 3 );
			fac_down = 3 * 441;
		}
#endif
		/*L_frame = ((L_frame_int*fac_down)+(*frac_mem))/fac_up;*/
		WorkLen = ( ( L_FRAME_PLUS * 2 * fac_down ) + ( st->right.decim_frac ) ) / fac_up;
	}
	else {
		WorkLen = L_frame;
	}
	/*---------------------------------------------------------------------*
  * Initialize pointers to speech vector.                               *
  *                                                                     *
  *                     20ms     20ms     20ms     20ms    >=20ms       *
  *             |----|--------|--------|--------|--------|--------|     *
  *           past sp   div1     div2     div3     div4    L_NEXT       *
  *             <--------  Total speech buffer (L_TOTAL_PLUS)  -------->     *
  *        old_speech                                                   *
  *                  <----- present frame (L_FRAME_PLUS) ----->              *
  *                  |        <------ new speech (L_FRAME_PLUS) ------->     *
  *                  |        |                                         *
  *                speech     |                                         *
  *                         new_speech                                  *
  *---------------------------------------------------------------------*/
	new_speech = old_speech + L_OLD_SPEECH_ST;
	speech = old_speech + L_TOTAL_ST - L_FRAME_PLUS - L_A_MAX - L_BSP;
	synth = old_synth + M;
	/*-----------------------------------------------------------------*
  * STEREO signal downsampling (codec working at 6.4kHz)            *
  * - decimate signal to fs=12.8kHz                                 *
  * - Perform 50Hz HP filtering of signal at fs=12.8kHz.            *
  *-----------------------------------------------------------------*/
	if( fscale == 0 ) {
		decim_12k8( channel_left, L_frame, sig_left, st->left.mem_decim, 0 );
		decim_12k8( channel_right, L_frame, sig_right, st->right.mem_decim, 0 );
		nb_samp = L_frame;
	}
	else {
		decim_split_12k8( channel_left, L_frame, sig_left, channel_left, L_FRAME_PLUS, fac_fs, fac_up, fac_down, WorkLen, st->left.mem_decim, &( st->left.decim_frac ) );

		nb_samp = decim_split_12k8( channel_right, L_frame, sig_right, channel_right, L_FRAME_PLUS, fac_fs, fac_up, fac_down, WorkLen, st->right.mem_decim, &( st->right.decim_frac ) );
	}
	hp50_12k8( sig_left, L_FRAME_PLUS, st->left.mem_sig_in, fscale );
	hp50_12k8( sig_right, L_FRAME_PLUS, st->right.mem_sig_in, fscale );
	/*-----------------------------------------------------------------*
  * Encode MONO low frequency band.                                 *
  * - Mix left and right channels (mono signal)                     *
  * - Perform fixed preemphasis through 1 - g z^-1                  *
  * - Encode low frequency band using ACELP/TCX model               *
  *-----------------------------------------------------------------*/
	mix_ch( sig_left, sig_right, new_speech, L_FRAME_PLUS, 1.0f, 1.0f );
	/* Apply preemphasis (for core codec only */
	E_UTIL_f_preemph( new_speech, PREEMPH_FAC, L_FRAME_PLUS, &( st->right.mem_preemph ) );
	/* copy memory into working space */
	mvr2r( st->old_speech_pe, old_speech, L_OLD_SPEECH_ST );
	mvr2r( st->old_synth, old_synth, M );
	if( use_case_mode == USE_CASE_B ) {
		for( i = 0; i < 4; i++ ) {
			st->stClass->vadFlag[i] = wb_vad( st->vadSt, &new_speech[256 * i] );
			if( st->stClass->vadFlag[i] == 0 ) {
				st->vad_hist++;
			}
			else {
				st->vad_hist = 0;
			}
			excType[i] = (short)( classifyExcitation( st->stClass, st->vadSt->level, (short)i ) );
		}
		/* encode mono lower band */
		coder_lf_b( codec_mode, speech, synth, mod, AqLF, st->window, param, ol_gain, T_out, p_out, excType, fscale, st );
	}
	else {
		for( i = 0; i < 4; i++ ) {
			excType[i] = 0;
		}
		/* encode mono lower band */
		coder_lf( codec_mode, speech, synth, mod, AqLF, st->window, param, ol_gain, T_out, p_out, excType, fscale, st );
	}
	for( i = 0; i < 4; i++ ) {
		mod[i] = excType[i];
	}
	/* update lower band memory for next frame */
	mvr2r( &old_speech[L_FRAME_PLUS], st->old_speech_pe, L_OLD_SPEECH_ST );
	mvr2r( &old_synth[L_FRAME_PLUS], st->old_synth, M );
	/*------------------------------------------------------------*
  * STEREO Bandwidth extension (2 channels used)               *
  * - fold and decimate higher band into new_speech_hf         *
  *   (2000Hz..6400Hz <-- 6400Hz..10800 Hz)                    *
  * - encode HF using 0.8kbps per channel.                     *
  *------------------------------------------------------------*/
	if( L_frame > L_FRAME8k ) {
		float  old_speech_hf[L_TOTAL_ST];
		float *new_speech_hf, *speech_hf;

		new_speech_hf = old_speech_hf + L_OLD_SPEECH_ST;
		speech_hf = old_speech_hf + L_TOTAL_ST - L_FRAME_PLUS - L_A_MAX - L_BSP;

		//left
		mvr2r( st->left.old_speech_hf, old_speech_hf, L_OLD_SPEECH_ST );
		if( fscale == 0 ) {
			decim_12k8( channel_left, L_frame, new_speech_hf, st->left.mem_decim_hf, ( fscale == 0 ) ? 1 : 2 );
		}
		else {
			/* left HF was stored into channel_left */
			mvr2r( channel_left, new_speech_hf, L_FRAME_PLUS );
		}
		mvr2r( &old_speech_hf[L_FRAME_PLUS], st->left.old_speech_hf, L_OLD_SPEECH_ST );
		if( StbrMode < 0 )
			mvr2r( speech_hf, channel_left, L_FRAME_PLUS + L_A_MAX + L_BSP );
		else
			coder_hf( mod, AqLF, speech, speech_hf, syn_hf, st->window, prm_hf_left, fscale, &( st->left ) );

		//right
		mvr2r( st->right.old_speech_hf, old_speech_hf, L_OLD_SPEECH_ST );
		if( fscale == 0 ) {
			decim_12k8( channel_right, L_frame, new_speech_hf, st->right.mem_decim_hf, ( fscale == 0 ) ? 1 : 2 );
		}
		else {
			/* right HF was stored into channel_right */
			mvr2r( channel_right, new_speech_hf, L_FRAME_PLUS );
		}
		if( StbrMode < 0 ) {
			for( i = 0; i < L_FRAME_PLUS + L_A_MAX + L_BSP; i++ )
				speech_hf[i] = 0.5f * speech_hf[i] + 0.5f * channel_left[i];
		}
		mvr2r( &old_speech_hf[L_FRAME_PLUS], st->right.old_speech_hf, L_OLD_SPEECH_ST );
		coder_hf( mod, AqLF, speech, speech_hf, syn_hf, st->window, prm_hf_right, fscale, &( st->right ) );
	}
	else {
		for( i = 0; i < NB_DIV * NPRM_BWE_DIV; i++ ) {
			prm_hf_right[i] = 0;
			prm_hf_left[i] = 0;
		}
	}
	/*------------------------------------------------------------*
  * STEREO low frequency band encoder                          *
  * - Mix left and right channels (mono signal)                *
  * - split mono and right signal (2k, hi).                    *
  * - perform parametric stereo encoding				   *
  *------------------------------------------------------------*/
	{
		float  old_mono_2k[L_TOTAL_ST_2k];
		float *new_mono_2k = old_mono_2k + L_OLD_SPEECH_2k;
		float *mono_2k = old_mono_2k + L_TOTAL_ST_2k - L_FRAME_2k - L_A_2k - L_FDEL_2k;
		float  old_chan_2k[L_TOTAL_ST_2k];
		float *new_chan_2k = old_chan_2k + L_OLD_SPEECH_2k;
		float *chan_2k = old_chan_2k + L_TOTAL_ST_2k - L_FRAME_2k - L_A_2k - L_FDEL_2k;
		float  old_mono_hi[L_TOTAL_ST_hi];
		float *new_mono_hi = old_mono_hi + L_OLD_SPEECH_hi;
		float *mono_hi = old_mono_hi + L_TOTAL_ST_hi - L_FRAME_PLUS - L_A_MAX;
		float  old_chan_hi[L_TOTAL_ST_hi];
		float *new_chan_hi = old_chan_hi + L_OLD_SPEECH_hi;
		float *chan_hi = old_chan_hi + L_TOTAL_ST_hi - L_FRAME_PLUS - L_A_MAX;
		/* copy memory into working space */
		mvr2r( st->old_speech, old_speech, L_OLD_SPEECH_ST );     /*528*/
		mvr2r( st->old_speech_2k, old_mono_2k, L_OLD_SPEECH_2k ); /*140*/
		mvr2r( st->old_speech_hi, old_mono_hi, L_OLD_SPEECH_hi ); /*448*/
		mix_ch( sig_left, sig_right, new_speech, L_FRAME_PLUS, 1.0f, 1.0f );
		/* do the lo,hi band-splitting on the mono signal */
		band_split_taligned_2k( new_speech, new_mono_2k, new_mono_hi, L_FRAME_PLUS ); /*128 samples delay  @ 12k8 */
		/* copy working space into memory */
		mvr2r( old_speech + L_FRAME_PLUS, st->old_speech, L_OLD_SPEECH_ST );
		mvr2r( old_mono_2k + L_FRAME_2k, st->old_speech_2k, L_OLD_SPEECH_2k );
		mvr2r( old_mono_hi + L_FRAME_PLUS, st->old_speech_hi, L_OLD_SPEECH_hi );
		/* copy memory into working space */
		mvr2r( st->old_chan, old_speech, L_OLD_SPEECH_ST );
		mvr2r( st->old_chan_2k, old_chan_2k, L_OLD_SPEECH_2k );
		mvr2r( st->old_chan_hi, old_chan_hi, L_OLD_SPEECH_hi );
		mvr2r( sig_right, new_speech, L_FRAME_PLUS );
		/* do the lo,hi band-splitting on the mono signal */
		band_split_taligned_2k( new_speech, new_chan_2k, new_chan_hi, L_FRAME_PLUS );
		/* copy working space into memory */
		mvr2r( old_speech + L_FRAME_PLUS, st->old_chan, L_OLD_SPEECH_ST );
		mvr2r( old_chan_2k + L_FRAME_2k, st->old_chan_2k, L_OLD_SPEECH_2k );
		mvr2r( old_chan_hi + L_FRAME_PLUS, st->old_chan_hi, L_OLD_SPEECH_hi );

		if( StbrMode < 0 ) {
			init_coder_stereo_x( st );
			st->mem_stereo_ovlp_size = L_OVLP_2k;
		}
		else
			coder_stereo_x( mono_hi, chan_hi, mono_2k, chan_2k, AqLF, StbrMode, prm_stereo, fscale, st );
	}
	/*--------------------------------------------------*
  * encode bits for serial stream                    *
  *--------------------------------------------------*/
	/* mode (0=ACELP 20ms, 1=TCX 20ms, 2=TCX 40ms, 3=TCX 80ms) */
	/* for 20-ms packetization, divide by 4 the 80-ms bitstream */
	nbits_pack = ( NBITS_CORE[codec_mode] + NBITS_BWE ) / 4;
	if( StbrMode >= 0 )
		nbits_pack += ( StereoNbits[StbrMode] + NBITS_BWE ) / 4;

	for( k = 0; k < NB_DIV; k++ ) {
		int2bin( mod[k], 2, &serial[k * nbits_pack] );
	}
	enc_prm( mod, codec_mode, param, serial, nbits_pack );
	if( StbrMode >= 0 ) {
		enc_prm_stereo_x( prm_stereo, serial, nbits_pack, NBITS_BWE, StbrMode );
		enc_prm_hf( mod, prm_hf_left, serial - NBITS_BWE / 4, nbits_pack );
	}
	enc_prm_hf( mod, prm_hf_right, serial, nbits_pack );

	return ( nb_samp );
}
/*-----------------------------------------------------------------*
 * Funtion  coder_amrwb_plus_mono                                  *
 *          ~~~~~~~~~~~~~~~~                                       *
 *   - Main mono coder routine.                                    *
 *                                                                 *
 *-----------------------------------------------------------------*/
int coder_amrwb_plus_mono(             /* output: number of sample processed */
    float             channel_right[], /* input: used on mono and stereo       */
    int               codec_mode,      /* input: AMR-WB+ mode (see cnst.h)         */
    int               L_frame,         /* input: 80ms frame size                   */
    short             serial[],        /* output: serial parameters                */
    Coder_State_Plus *st,              /* i/o : coder memory state                 */
    short             use_case_mode,
    int               fscale )
{
	/* LPC coefficients of lower frequency */
	float AqLF[( NB_SUBFR + 1 ) * ( M + 1 )];
	int   param[NB_DIV * NPRM_DIV];
	int   prm_hf_right[NB_DIV * NPRM_BWE_DIV];
	/* vector working at fs=12.8kHz */
	float  old_speech[L_TOTAL_ST];
	float  old_synth[M + L_FRAME_PLUS];
	float *speech, *new_speech;
	float *synth;
	/* Scalars */
	int i, k, nbits_pack;
	int mod[NB_DIV];
	/* LTP parameters for high band */
	float ol_gain[NB_DIV];
	int   T_out[NB_DIV];
	float p_out[NB_DIV];
	int   nb_samp, fac_fs;
	/*ClassB parameters*/
	short excType[4];
	int   WorkLen, fac_up, fac_down;

	if( ( L_frame - L_FRAME32k ) == 0 ) {
		fac_fs = FSCALE_DENOM * 3 / 2;
	}
	else {
		fac_fs = fscale;
	}
	/* 48k setting*/
	fac_up = ( fac_fs << 3 );
	fac_down = 180 * 8;

	if( fscale != 0 ) {

#ifdef FILTER_44kHz
		if( ( L_frame - ( 2 * L_FRAME44k ) ) == 0 ) {
			fac_up = ( fac_fs << 3 );
			fac_down = 3 * 441;
		}
#endif
		/*L_frame = ((L_frame_int*fac_down)+(*frac_mem))/fac_up;*/
		WorkLen = ( ( L_FRAME_PLUS * 2 * fac_down ) + ( st->right.decim_frac ) ) / fac_up;
	}
	else {
		WorkLen = L_frame;
	}
	/*---------------------------------------------------------------------*
  * Initialize pointers to speech vector.                               *
  *                                                                     *
  *                     20ms     20ms     20ms     20ms    >=20ms       *
  *             |----|--------|--------|--------|--------|--------|     *
  *           past sp   div1     div2     div3     div4    L_NEXT       *
  *             <--------  Total speech buffer (L_TOTAL_PLUS)  -------->     *
  *        old_speech                                                   *
  *                  <----- present frame (L_FRAME_PLUS) ----->              *
  *                  |        <------ new speech (L_FRAME_PLUS) ------->     *
  *                  |        |                                         *
  *                speech     |                                         *
  *                         new_speech                                  *
  *---------------------------------------------------------------------*/
	new_speech = old_speech + L_OLD_SPEECH_ST;
	speech = old_speech + L_TOTAL_ST - L_FRAME_PLUS - L_A_MAX - L_BSP;
	synth = old_synth + M;
	/*-----------------------------------------------------------------*
  * MONO/STEREO signal downsampling (codec working at 6.4kHz)       *
  * - decimate signal to fs=12.8kHz                                 *
  * - Perform 50Hz HP filtering of signal at fs=12.8kHz.            *
  * - Perform fixed preemphasis through 1 - g z^-1                  *
  * - Mix left and right channels into sum and difference signals   *
  * - perform parametric stereo encoding                            *
  *-----------------------------------------------------------------*/
	if( fscale == 0 ) {
		decim_12k8( channel_right, L_frame, new_speech, st->right.mem_decim, 0 );
		nb_samp = L_frame;
	}
	else {

		nb_samp = decim_split_12k8( channel_right, L_frame, new_speech, channel_right, L_FRAME_PLUS, fac_fs, fac_up, fac_down, WorkLen, st->right.mem_decim, &( st->right.decim_frac ) );
	}
	hp50_12k8( new_speech, L_FRAME_PLUS, st->right.mem_sig_in, fscale );
	/*------------------------------------------------------------*
  * Encode MONO low frequency band using ACELP/TCX model       *
  *------------------------------------------------------------*/
	/* Apply preemphasis (for core codec only */
	E_UTIL_f_preemph( new_speech, PREEMPH_FAC, L_FRAME_PLUS, &( st->right.mem_preemph ) );
	/* copy memory into working space */
	mvr2r( st->old_speech_pe, old_speech, L_OLD_SPEECH_ST );
	mvr2r( st->old_synth, old_synth, M );
	if( use_case_mode == USE_CASE_B ) {
		for( i = 0; i < 4; i++ ) {
			st->stClass->vadFlag[i] = wb_vad( st->vadSt, &new_speech[256 * i] );
			if( st->stClass->vadFlag[i] == 0 ) {
				st->vad_hist++;
			}
			else {
				st->vad_hist = 0;
			}
			excType[i] = (short)( classifyExcitation( st->stClass, st->vadSt->level, (short)i ) );
		}
		/* encode mono lower band */
		coder_lf_b( codec_mode, speech, synth, mod, AqLF, st->window, param, ol_gain, T_out, p_out, excType, fscale, st );
	}
	else {
		for( i = 0; i < 4; i++ ) {
			excType[i] = 0;
		}
		/* encode mono lower band */
		coder_lf( codec_mode, speech, synth, mod, AqLF, st->window, param, ol_gain, T_out, p_out, excType, fscale, st );
	}
	for( i = 0; i < 4; i++ ) {
		mod[i] = excType[i];
	}
	st->prev_mod = (short)mod[3];
	/* update lower band memory for next frame */
	mvr2r( &old_speech[L_FRAME_PLUS], st->old_speech_pe, L_OLD_SPEECH_ST );
	mvr2r( &old_synth[L_FRAME_PLUS], st->old_synth, M );
	/*------------------------------------------------------------*
  * MONO/STEREO Bandwidth extension (2 channels used in stereo)*
  * - fold and decimate higher band into new_speech_hf         *
  *   (2000Hz..6400Hz <-- 6400Hz..10800 Hz)                    *
  * - encode HF using 0.8kbps per channel.                     *
  *------------------------------------------------------------*/
	if( L_frame > L_FRAME8k ) {
		float  old_speech_hf[L_TOTAL_ST];
		float *new_speech_hf, *speech_hf;
		float *syn_hf = old_synth;
		new_speech_hf = old_speech_hf + L_OLD_SPEECH_ST;
		speech_hf = old_speech_hf + L_TOTAL_ST - L_FRAME_PLUS - L_A_MAX - L_BSP;
		mvr2r( st->right.old_speech_hf, old_speech_hf, L_OLD_SPEECH_ST );
		if( fscale == 0 ) {
			decim_12k8( channel_right, L_frame, new_speech_hf, st->right.mem_decim_hf, ( fscale == 0 ) ? 1 : 2 );
		}
		else {
			/* HF was stored into channel_right */
			mvr2r( channel_right, new_speech_hf, L_FRAME_PLUS );
		}
		mvr2r( &old_speech_hf[L_FRAME_PLUS], st->right.old_speech_hf, L_OLD_SPEECH_ST );
		coder_hf( mod, AqLF, speech, speech_hf, syn_hf, st->window, prm_hf_right, fscale, &( st->right ) );
	}
	else {
		for( k = 0; k < NB_DIV * NPRM_BWE_DIV; k++ ) {
			prm_hf_right[k] = 0;
		}
	}
	/*--------------------------------------------------*
  * encode bits for serial stream                    *
  *--------------------------------------------------*/
	/* mode (0=ACELP 20ms, 1=TCX 20ms, 2=TCX 40ms, 3=TCX 80ms) */
	/* for 20-ms packetization, divide by 4 the 80-ms bitstream */
	nbits_pack = ( NBITS_CORE[codec_mode] + NBITS_BWE ) / 4;
	for( k = 0; k < NB_DIV; k++ ) {
		int2bin( mod[k], 2, &serial[k * nbits_pack] );
	}
	enc_prm( mod, codec_mode, param, serial, nbits_pack );
	enc_prm_hf( mod, prm_hf_right, serial, nbits_pack );
	return ( nb_samp );
}
