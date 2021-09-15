#ifndef proto_func_h
#define proto_func_h

#include <stdio.h>

#include "../include/amr_plus.h"
#include "../include/mem.h"
#include "../lib_amr/typedef.h"

/* AMR-WB+ Function prototypes */

/* Decoder files */

/* tcx_ecu.c */
void adapt_low_freq_deemph_ecu( float xri[], int lg, Decoder_State_Plus *st );
void reconst_spect( float xri[], float old_xri[], int n_pack, int bfi[], int lg, int last_mode, float buf[] );

/* Encoder files */

/* lag window lag_wind.c */
void init_lag_wind( float bwe,    /* input : bandwidth expansion */
    float                 f_samp, /* input : sampling frequency */
    float                 wnc,    /* input : white noise correction factor */
    int                   m );                      /* input : order of LP filter */
void lag_wind( float r[],         /* in/out: autocorrelations */
    int              m );                      /* input : order of LP filter */

/* Common files */

/* Adaptive Low freq emphasis in alf_emph.c */
void adap_low_freq_emph( float xri[], int lg );
void adap_low_freq_deemph( float xri[], int lg );

/* bit packing and unpacking functions in bits.c */
int  bin2int(                 /* output: recovered integer value */
    int    no_of_bits,       /* input : number of bits associated with value */
    short *bitstream );      /* input : address where bits are read */
void int2bin( int value,      /* input : value to be converted to binary */
    int           no_of_bits, /* input : number of bits associated with value */
    short *       bitstream );       /* output: address where bits are written */

int over_fs(         /* number of sample oversampled       */
    float sig_in[],  /* input:  signal to oversample       */
    float sig_out[], /* output: signal oversampled         */
    int   lg,        /* input:  length of output           */
    int   fac_down,  /* input:  fs*12/fac_down = 44k/48k   */
    float mem[],     /* in/out: mem[2*L_FILT_OVER_FS]      */
    int * frac_mem   /* in/out: interpol fraction memory   */
);
int decim_fs(        /* number of sample decimated         */
    float sig_in[],  /* input:  signal to decimate         */
    int   lg,        /* input:  length of input            */
    float sig_out[], /* output: signal decimated           */
    int   fac_up,    /* input:  44k/48k *fac_up/12 = fs    */
    float mem[],     /* in/out: mem[2*L_FILT_DECIM_FS]     */
    int * frac_mem   /* in/out: interpol fraction memory   */
);
/* Resampling routines in deci12k8.c */
int decim_split_12k8(   /* number of sample decimated         */
    float sig_fs[],     /* input:  signal to decimate         */
    int   lg_input,     /* input:  2*L_FRAME44k if 44kHz      */
    float sig12k8_lf[], /* output: LF decimated signal        */
    float sig12k8_hf[], /* output: HF decimated signal        */
    int   lg,           /* input:  length of LF and HF        */
    int   fac_fs,       /* input:  >=32 (32 = base fs)        */
    int   fac_up,       /* (i)  :  */
    int   fac_down,     /* (i)  :  */
    int   L_frame,      /* (i)  :  */
    float mem[],        /* in/out: mem[L_MEM_DECIM_SPLIT]     */
    int * frac_mem );    /* in/out: interpol fraction memory   */

int join_over_12k8(     /* number of sample oversampled       */
    float sig12k8_lf[], /* input:  LF signal (fs=12k8)        */
    float sig12k8_hf[], /* input:  HF signal (fs=12k8)        */
    int   lg,           /* input:  length of LF and HF        */
    float sig_fs[],     /* output: oversampled signal         */
    int   lg_output,    /* input:  L_FRAME44k if 44kHz        */
    int   fac_fs,       /* input:  >=32 (32 = base fs)        */
    float mem[],        /* in/out: mem[L_MEM_JOIN_OVER]       */
    int * frac_mem );    /* in/out: interpol fraction memory   */

void decim_12k8( float sig_fs[],     /* input: signal to decimate */
    int                lg,           /* input: length of input */
    float              sig12k8[],    /* output: decimated signal */
    float              mem[],        /* in/out: memory (2*L_FILT_FS) */
    int                band );                      /* input: 0=0..6.4k, 1=6.4..10.8k */
void oversamp_12k8( float sig12k8[], /* input: signal to oversampling */
    float                 sig_fs[],  /* output: oversampled signal */
    int                   lg,        /* input: length of output */
    float                 mem[],     /* in/out: memory (2*L_FILT) */
    int                   band,      /* input: 0=0..6.4k, 1=6.4..10.8k */
    int                   add );
void interpol( float *signal, float *signal_int, int L_frame_int, const float *filter, int nb_coef, int fac_up, int fac_down, float gain );

/* FFT routines in fft3.c and fft9.c */
void fft3( float X[], float Y[], short n );
void ifft3( float X[], float Y[], short n );
void fft9( float X[], float Y[], short n );
void ifft9( float Y[], float X[], short n );
void fft_rel( float x[], short n, short m );
void ifft_rel( float x[], short n, short m );

/* TCX gain functions in gaintcx.c */
int   q_gain_tcx(     /* output: return quantization index */
    float  xnq[],   /* (i) : quantized vector */
    int    lg,      /* (i) : frame size */
    float *gain );  /* in/out: quantized gain */
float d_gain_tcx(     /* output: gain */
    int    index,     /* (i) : index */
    float  code[],    /*(i) : quantized vector */
    int    lcode,     /* (i) : frame size */
    int    bfi,       /* (i) : 1=gain lost */
    float *old_rms ); /* (i/o): for frame recovery */

/* hf_func.c */
float match_gain_6k4( float *AqLF, float *AqHF );
void  int_gain( float old_gain, float new_gain, float *gain, int nb_subfr );
void  soft_exc_hf( float *exc_hf, float *mem );
void  soft_exc_hf_new( float *exc_hf, float *mem, int l_frame );
void  smooth_ener_hf( float *HF, float *threshold );

/* High-pass filter in hp50.c */
void hp50_12k8( Float32 signal[], Word32 lg, Float32 mem[], Word32 fscale );

/* int_lpc.c */
void int_lpc_np1( float isf_old[], /* input : LSFs from past frame */
    float               isf_new[], /* input : LSFs from present frame */
    float               a[],       /* output: LP coefficients in both subframes */
    int                 nb_subfr,  /* input: number of subframe */
    int                 m );                       /* input : order of LP filter */

/* pitch predictor in pit_fr4.c */
void pred_lt4( float exc[], /* in/out: excitation buffer */
    int              T0,    /* input : integer pitch lag */
    int              frac,  /* input : fraction of lag */
    int              L_subfr );          /* input : subframe size */

/* q_gn_hf.c */
void q_gain_hf( float *gain,   /* input : gain of 4 subfr */
    float *            gain_q, /* output: quantized gains */
    int *              indice );             /* output: indices */

void d_gain_hf( int indice, /* input: quantization indices */
    float *         gain_q, /* output: quantized gains */
    float *         past_q, /* i/o : past quantized gain (1 word) */
    int             bfi );              /* input : Bad frame indicator */

/* q_isf_hf.c */
void q_isf_hf( float *isf1,   /* input : ISF in the frequency domain (0..6400) */
    float *           isf_q,  /* output: quantized ISF */
    float *           past_q, /* i/o : past quantized isf (for MA prediction) */
    int *             indice, /* output: quantization indices (7 words) */
    const float *     mean_isf_hf,
    const float *     dico1_isf_hf );

void d_isf_hf( int *indice, /* input: quantization indices */
    float *         isf_q,  /* output: quantized ISFs in the cosine domain */
    float *         past_q, /* i/o : past quantized isf (for MA prediction) */
    int             bfi,    /* input : Bad frame indicator */
    const float *   mean_isf_hf,
    const float *   dico1_isf_hf );

/* qpisf_2s.c */
void qpisf_2s_46b( float *isf1,      /* input : ISF in the frequency domain (0..6400) */
    float *               isf_q,     /* output: quantized ISF */
    float *               past_isfq, /* i/0 : past ISF quantizer */
    int *                 indice,    /* output: quantization indices (7 words) */
    int                   nb_surv );                   /* input : number of survivor (1, 2, 3 or 4) */
void dpisf_2s_46b( int *indice,      /* input: quantization indices */
    float *             isf_q,       /* output: quantized ISFs in the cosine domain */
    float *             past_isfq,   /* i/0 : past ISF quantizer */
    float *             isfold,      /* input : past quantized ISF */
    float *             isf_buf,     /* input : isf buffer */
    int                 bfi,         /* input : Bad frame indicator */
    int                 bfi_2nd_st,  /* input : 2nd stage bfi mask (bin: 011111) */
    int                 enc_dec );

/* RE8 lattice quantiser functions in re8_*.c */
void RE8_PPV( float x[], int y[] );
void RE8_cod( int *y, int *n, long *I, int *k ); /* encoder only */
void RE8_dec( int nq, long I, int kv[], int y[] );
void RE8_vor( int y[], int *n, int k[], int c[], int *ka );
void re8_coord( int y[], int k[] );
void re8_k2y( int k[], int m, int y[] );

/* rnd_ph16.c */
void rnd_ph16( short *seed, float *xri, int lg );

/* util.c */
void set_zero( float *x, int n );
void mvr2r( float x[], /* input : input vector */
    float         y[], /* output: output vector */
    int           n );           /* input : vector size */
void mvs2s( short x[], /* input : input vector */
    short         y[], /* output: output vector */
    int           n );           /* input : vector size */
void mvi2i( int x[],   /* input : input vector */
    int         y[],   /* output: output vector */
    int         n );           /* input : vector size */
void mvr2s( float x[], /* input : input vector */
    short         y[], /* output: output vector */
    int           n );           /* input : vector size */
void mvs2r( short x[], /* input : input vector */
    float         y[], /* output: output vector */
    int           n );           /* input : vector size */
void pessimize();

/* util_stereo_x.c */
#include "util_stereo_x.h"
float my_max( float x, float y );
float my_min( float x, float y );
void  pmsvq( float *y, int **prm, float *x, float *old_x, const PMSVQ *filt_hi_pmsvq );
void  pmsvq_inv( float *y, int **prm, float *old_y, const int bfi, const PMSVQ *filt_hi_pmsvq );
void  pvq( float *x, float *old_x, float *cb, float *cbm, float a, int n, int cb_size, int *winner );
void  syn_filt( float a[],   /* input : LP filter coefficients */
     int              m,     /* input : order of LP filter */
     float            x[],   /* input : input signal */
     float            y[],   /* output: output signal */
     int              l,     /* input : size of filtering */
     float            mem[], /* in/out: initial filter states */
     int              update_m );         /* input : update memory flag: 0 --> no memory update 1 --> update */
void  residu( float *a,      /* input : LP filter coefficients */
     int             m,      /* input : order of LP filter */
     float *         x,      /* input : input signal (usually speech) */
     float *         y,      /* output: output signal (usually residual) */
     int             l );                /* input : size of filtering */
int   cholsolc( float r[HI_FILT_ORDER][HI_FILT_ORDER], float c[HI_FILT_ORDER], float h[HI_FILT_ORDER], int n );
float glev_s( float *b,         /* output: filter coefficients */
    float *          r,         /* input : vector of autocorrelations */
    float *          z,         /* input: vector of cross correlations */
    int              m );                    /* input : order of LP filte */
void  crosscorr( float *vec1,   /* (i) : Input vector 1 */
     float *            vec2,   /* (i) : Input vector 2 */
     float *            result, /* (o) : Output result vector */
     int                length, /* (i) : Length of input vectors */
     int                minlag, /* (i) : Minimum lag */
     int                maxlag );              /* (i) : Maximum lag */

/* window.c */

/* File IO interface routines: from
 * read_dat.c, wavefiletools.c and writ_dat.c */
int   read_data(               /* return: number of data successfully read */
    FILE *fp,                /* input : data file (16-bit words) */
    float data[],            /* output: speech data */
    int   size );              /* input : number of samples */
void  writ_data( float data[], /* input : data */
     int               size,   /* input : number of samples */
     FILE *            fp );               /* output: file pointer */
FILE *Wave_fopen( char *Filename, char *Mode, short *NumOfChannels, long *SamplingRate, short *BitsPerSample, long *DataSize );
void  Wave_fclose( FILE *FilePtr, short BitsPerSample );

/*---------------------------------------------------------------------*
 *              main routines                                          *
 *---------------------------------------------------------------------*/

void init_coder_amrwb_plus( Coder_State_Plus *st, int num_chan, int fscale, short use_case_mode, short full_reset );

int coder_amrwb_plus_stereo( float channel_right[], /* input: used on mono and stereo */
    float                          channel_left[],  /* input: used on stereo only */
    int                            codec_mode,      /* input: AMR-WB+ mode (see cnst.h) */
    int                            L_frame,         /* input: 80ms frame size */
    short                          serial[],        /* output: serial parameters */
    Coder_State_Plus *             st,              /* i/o : coder memory state */
    short                          useCaseB,
    int                            bwe_flag, /* 32kHz NBWE */
    int                            StbrMode );

int coder_amrwb_plus_mono( float channel_right[], /* input: used on mono and stereo */
    int                          codec_mode,      /* input: AMR-WB+ mode (see cnst.h) */
    int                          L_frame,         /* input: 80ms frame size */
    short                        serial[],        /* output: serial parameters */
    Coder_State_Plus *           st,              /* i/o : coder memory state */
    short                        useCaseB,
    int                          bwe_flag /* 32kHz NBWE */
);

void coder_amrwb_plus_mono_first( float channel_right[], /* input: used on mono and stereo */
    int                                 n_channel,       /* input: 1 or 2 (mono/stereo) */
    int                                 L_frame,         /* input: frame size */
    int                                 L_next,          /* input: lookahead */
    int                                 bwe_flag,        /* for 32kHz NBWE */
    Coder_State_Plus *                  st               /* i/o : coder memory state */
);

int coder_amrwb_plus_first( float channel_right[], /* input: used on mono and stereo */
    float                         channel_left[],  /* input: used on stereo only */
    int                           n_channel,       /* input: 1 or 2 (mono/stereo) */
    int                           L_frame,         /* input: frame size */
    int                           L_next,          /* input: lookahead */
    int                           bwe_flag,        /* AriL: for 32kHz NBWE */
    Coder_State_Plus *            st               /* i/o : coder memory state */
);

void init_decoder_amrwb_plus( Decoder_State_Plus *st, int num_chan, int fscale, short full_reset );

int decoder_amrwb_plus( int codec_mode,      /* input: AMR-WB+ mode (see cnst.h) */
    short                   serial[],        /* input: serial parameters (4x20ms) */
    int                     bad_frame[],     /* input: bfi (bad_frame[4]) */
    int                     L_frame,         /* input: frame size of synthesis */
    int                     n_channel,       /* input: 1 or 2 (mono/stereo) */
    float                   channel_right[], /* (o): used on mono and stereo */
    float                   channel_left[],  /* (o): used on stereo only */
    Decoder_State_Plus *    st,              /* i/o : decoder memory state */
    int                     fscale,
    int                     StbrMode,
    int                     mono_dec_stereo,
    short                   upsamp_fscale );

void decoder_amrwb_plus_1( float *chan_right,
    float *                       chan_left,
    int *                         mod,
    int *                         param,
    int *                         prm_hf_right,
    int *                         prm_hf_left,
    int *                         nbits_AVQ,
    int                           codec_mode,
    int *                         bad_frame,
    int *                         bad_frame_hf,
    float *                       AqLF,
    float *                       synth,
    int *                         pitch,
    float *                       pit_gain,
    Decoder_State_Plus *          st,
    int                           n_channel,
    int                           L_frame,
    int                           bwe_flag,
    int                           mono_dec_stereo );

/*---------------------------------------------------------------------*
 *              low freq band routines (0..6400Hz)                     *
 *---------------------------------------------------------------------*/

void init_coder_lf( Coder_State_Plus *st );
void coder_lf( int    codec_mode,   /* (i) : AMR-WB+ mode (see cnst.h) */
    float             speech[],     /* (i) : speech vector [-M..L_FRAME_PLUS+L_NEXT] */
    float             synth[],      /* (o) : synthesis vector [-M..L_FRAME_PLUS] */
    int               mod[],        /* (o) : mode for each 20ms frame (mode[4] */
    float             AqLF[],       /* (o) : quantized coefficients (AdLF[16]) */
    float             window[],     /* (i) : window for LPC analysis */
    int               param[],      /* (o) : parameters (NB_DIV*NPRM_DIV) */
    float             ol_gain[],    /* (o) : open-loop LTP gain */
    int               ave_T_out[],  /* (o) : average LTP lag */
    float             ave_p_out[],  /* (o) : average LTP gain */
    short             coding_mod[], /* (i) : selected mode for each 20ms */
    int               fscale,
    Coder_State_Plus *st /* i/o : coder memory state */
);
void coder_lf_b( int  codec_mode,   /* (i) : AMR-WB+ mode (see cnst.h) */
    float             speech[],     /* (i) : speech vector [-M..L_FRAME_PLUS+L_NEXT] */
    float             synth[],      /* (o) : synthesis vector [-M..L_FRAME_PLUS] */
    int               mod[],        /* (o) : mode for each 20ms frame (mode[4] */
    float             AqLF[],       /* (o) : quantized coefficients (AdLF[16]) */
    float             window[],     /* (i) : window for LPC analysis */
    int               param[],      /* (o) : parameters (NB_DIV*NPRM_DIV) */
    float             ol_gain[],    /* (o) : open-loop LTP gain */
    int               ave_T_out[],  /* (o) : average LTP lag */
    float             ave_p_out[],  /* (o) : average LTP gain */
    short             coding_mod[], /* (i) : selected mode for each 20ms */
    int               pit_adj,
    Coder_State_Plus *st /* i/o : coder memory state */
);

void init_decoder_lf( Decoder_State_Plus *st );

void decoder_lf( int    mod[],       /* (i) : mode for each 20ms frame (mode[4] */
    int                 prm[],       /* (i) : parameters */
    int                 nbits_AVQ[], /* (i) : for each frame (nbits_AVQ[4]) */
    int                 codec_mode,  /* (i) : AMR-WB+ mode (see cnst.h) */
    int                 bad_frame[], /* (i) : for each frame (bad_frame[4]) */
    float               AqLF[],      /* (o) : decoded coefficients (AdLF[16]) */
    float               fexc[],      /* (o) : decoded excitation */
    float               fsynth[],    /* (o) : decoded synthesis */
    int                 pitch[],     /* (o) : decoded pitch (pitch[16]) */
    float               pit_gain[],  /* (o) : decoded pitch gain (pit_gain[16]) */
    int                 fscale,
    Decoder_State_Plus *st ); /* i/o : coder memory state */

void init_bass_postfilter( Decoder_State_Plus *st );
void bass_postfilter( float *synth_in,  /* (i) : 12.8kHz synthesis to postfilter */
    int *                    T_sf,      /* (i) : Pitch period for all subframe (T_sf[16]) */
    float *                  gainT_sf,  /* (i) : Pitch gain for all subframe (gainT_sf[16]) */
    float *                  synth_out, /* (o) : filtered synthesis (with delay=L_SUBFR+L_FILT) */
    int                      fscale,
    Decoder_State_Plus *     st ); /* i/o : decoder memory state */

/*---------------------------------------------------------------------*
 *              high freq band routines (6400Hz..10800Hz)              *
 *---------------------------------------------------------------------*/

void init_coder_hf( Coder_StState *st );
void coder_hf( int mod[],       /* (i) : mode for each 20ms frame (mode[4] */
    float          AqLF[],      /* (i) : Q coeff in lower band (AdLF[16]) */
    float          speech[],    /* (i) : speech vector [-M..L_FRAME_PLUS] */
    float          speech_hf[], /* (i) : HF speech vec [-MHF..L_FRAME_PLUS+L_NEXT] */
    float          synth_hf[],  /* (o) : HF synthesis [0..L_FRAME_PLUS] */
    float          window[],    /* (i) : window for LPC analysis */
    int            param[],     /* (o) : parameters (NB_DIV*NPRM_HF_DIV) */
    int            fscale,
    Coder_StState *st ); /* i/o : coder memory state */

void init_decoder_hf( Decoder_StState *st );
void decoder_hf( int mod[],           /* (i) : mode for each 20ms frame (mode[4] */
    int              param[],         /* (i) : parameters */
    int              param_other[],   /* (i) : parameters for the right channel in case of mono decoding with stereo bitstream */
    int              mono_dec_stereo, /* 1=Mono decoding with stereo bitstream */
    int              bad_frame[],     /* (i) : for each frame (bad_frame[4]) */
    float            AqLF[],          /* (i) : decoded coefficients (AdLF[16]) */
    float            exc[],           /* (i) : decoded excitation */
    float            fsynth_hf[],     /* (o) : decoded synthesis */
    float            mem_lpc_hf[],
    float *          mem_gain_hf,
    short *          ramp_state,
    Decoder_StState *st /* i/o : decoder memory state */
);

/*---------------------------------------------------------------------*
 *                   stereo routines (0..6400Hz)                       *
 *---------------------------------------------------------------------*/

/*-----------------------------------------------------------------*
 * Function coder_stereo                                            *
 * ~~~~~~~~~~~~~~~~~~~~                                            *
 *   ->Principle stereo coder routine (working at fs=12.8kHz).     *
 *                                                                 *
 * Note: HF band are encoded twice (2 channels) using 0.8kbps BWE. *
 *       Usage of 2xBWE for stereo provide better time domain      *
 *       stereo definition in HF without increasing the bit-rate.  *
 *       Another advantage is that the stereo decoder is limited   *
 *       to the lower band (fs=12.8kHz) and this reduce the        *
 *       overall complexity of the AMR-WB+ codec.  Also, this      *
 *       solution is not dependent of the AMR-WB+ mode where many  *
 *       different sampling frequencies are used (16, 24, 32 kHz). *
 *-----------------------------------------------------------------*/

void coder_stereo_x( float speech_hi[], /* (i) : Mixed channel, hi */
    float                  chan_hi[],   /* (i) : Right channel, hi */
    float                  speech_2k[], /* (i) : Mixed channel, lo */
    float                  chan_2k[],   /* (i) : Right channel, lo */
    float                  AqLF[],
    int                    brMode,
    int                    param[], /* (o) : Encoded parameters */
    int                    fscale,
    Coder_State_Plus *     st ); /* (i/o): Encoder states */

/*-----------------------------------------------------------------*
 * Function  init_coder_stereo                                      *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~                                      *
 *   ->Initialization of variables for the stereo coder.           *
 *-----------------------------------------------------------------*/
void init_coder_stereo_x( Coder_State_Plus *st ); /* (i/o) : Encoder states. */

/*-----------------------------------------------------------------*
 * Function  end_prm_stereo									       *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~                                      *
 *   ->Encoding of stereo parameters.					           *
 *-----------------------------------------------------------------*/

void enc_prm_stereo( int param[],    /* (i) : parameters */
    short                serial[],   /* (o) : serial bits stream */
    int                  nbits_pack, /* (i) : number of bits per packet of 20ms */
    int                  brMode,
    Coder_State_Plus *   st /* (i) : Encoder states */
);

void enc_prm_stereo_x( int param[],    /* (i) : parameters */
    short                  serial[],   /* (o) : serial bits stream */
    int                    nbits_pack, /* (i) : number of bits per packet of 20ms */
    int                    nbits_bwe,  /* (i) : number of BWE bits per 20ms  */
    int                    brMode );

void filt_design( float *spL,    /* (i) : Speech input (L) */
    float *              spR,    /* (i) : Speech input (R) */
    float *              spM,    /* (i) : Speech input (M) */
    float *              filter, /* (o) : Filter coeff output */
    Coder_State_Plus *   st,     /* (i/o): Filter states */
    int                  frame_length,
    int                  lb,             /* (i) : Look-back */
    int                  anaframe_length /* (i) : Analysis frame length */
);

void enc_quant_frame( int *prm, short *ptr, Coder_State_Plus *st );

void init_decoder_stereo_x( Decoder_State_Plus *st );

void decoder_stereo_x( int param[],
    int                    bad_frame[],
    float                  sig_left[],
    float                  sig_right[],
    float                  AqLF[],
    int                    StbrMode,
    int                    fscale,
    Decoder_State_Plus *   st );

void conceal_filt( float *filt_lo_ptr, float *filt_hi_ptr, Decoder_State_Plus *st );

void dec_quant_frame( int *prm, short *ptr, Decoder_State_Plus *st );

/*---------------------------------------------------------------------*
 *             Parameters encoding/decoding routines                   *
 *---------------------------------------------------------------------*/

void enc_prm( int mode[],     /* (i) : frame mode (mode[4], 4 division) */
    int           codec_mode, /* (i) : AMR-WB+ mode (see cnst.h) */
    int           param[],    /* (i) : parameters */
    short         serial[],   /* (o) : serial bits stream */
    int           nbits_pack  /* (i) : number of bits per packet of 20ms */
);
void enc_prm_hf( int mod[],     /* (i) : frame mode (mode[4], 4 division) */
    int              param[],   /* (i) : parameters */
    short            serial[],  /* (o) : serial bits stream */
    int              nbits_pack /* (i) : number of bits per packet of 20ms */
);
void enc_prm_stereo( int param[],    /* (i) : parameters */
    short                serial[],   /* (o) : serial bits stream */
    int                  nbits_pack, /* (i) : number of bits per packet of 20ms */
    int                  brMode,
    Coder_State_Plus *   st /* (i) : Encoder states */
);

void enc_prm_stereo_x( int param[],    /* (i) : parameters */
    short                  serial[],   /* (o) : serial bits stream */
    int                    nbits_pack, /* (i) : number of bits per packet of 20ms */
    int                    nbits_bwe,  /* (i) : number of BWE bits per 20ms  */
    int                    brMode );

void dec_prm( int mod[],       /* (i) : frame mode (mode[4], 4 frames) */
    int           bad_frame[], /* (i) : bfi for 4 frames (bad_frame[4]) */
    short         serial[],    /* (i) : serial bits stream */
    int           nbits_pack,  /* (i) : number of bits per packet of 20ms */
    int           codec_mode,  /* (i) : AMR-WB+ mode (see cnst.h) */
    int           param[],     /* (o) : decoded parameters */
    int           nbits_AVQ[] );         /* (o) : nb of bits for AVQ (4 division) */

void dec_prm_hf( int mod[],       /* (i) : frame mode (mode[4], 4 frames) */
    int              bad_frame[], /* (i) : bfi for 4 frames (bad_frame[4]) */
    short            serial[],    /* (i) : serial bits stream */
    int              nbits_pack,  /* (i) : number of bits per packet of 20ms */
    int              param[] );                /* (o) : decoded parameters */

void dec_prm_stereo( int bad_frame[], /* (i) : bfi for 4 frames (bad_frame[4]) */
    short                serial[],    /* (i) : serial bits stream */
    int                  nbits_pack,  /* (i) : number of bits per packet of 20ms */
    int                  param[],     /* (o) : decoded parameters */
    int                  brMode,
    Decoder_State_Plus * st );

void dec_prm_stereo_x( int bad_frame[], /* (i) : bfi for 4 frames (bad_frame[4]) */
    short                  serial[],    /* (i) : serial bits stream */
    int                    nbits_pack,  /* (i) : number of bits per packet of 20ms */
    int                    nbits_bwe,   /* (i) : number of BWE bits per 20ms  */
    int                    param[],     /* (o) : decoded parameters */
    int                    brMode,
    Decoder_State_Plus *   st );

/*---------------------------------------------------------------------*
 *              ACELP routines                                         *
 *---------------------------------------------------------------------*/

void coder_acelp( float A[],        /* input: coefficients 4xAz[M+1] */
    float               Aq[],       /* input: coefficients 4xAz_q[M+1] */
    float               speech[],   /* input: speech[-M..lg] */
    float *             mem_wsp,    /* in/out: wsp memory */
    float *             mem_wsyn,   /* in/out: wsyn memory */
    float               synth[],    /* in/out: synth[-M..lg] */
    float               exc[],      /* in/out: exc[-(PIT_MAX+L_INTERPOL)..lg+1] */
    float               wovlp[],    /* in/out: wovlp[0..lg+128] */
    int                 lg,         /* input: frame length */
    int                 codec_mode, /* input: AMR_WB+ mode (see cnst.h) */
    float               norm_corr,
    float               norm_corr2,
    int                 T_op,    /* input: open-loop LTP */
    int                 T_op2,   /* input: open-loop LTP */
    int                 T_out[], /* output: integer pitch-lag */
    float               p_out[], /* output: pitch gain */
    float               c_out[], /* output: fixed codebook gain */
    int                 fscale,
    int *               prm ); /* output: acelp parameters */

void decoder_acelp(
    int                 prm[],      /* input: parameters               */
    float               A[],        /* input: coefficients NxAz[M+1]   */
    int                 lg,         /* input: frame length             */
    int                 codec_mode, /* input: AMR-WB+ mode (see cnst.h)*/
    int                 bfi,        /* input: 1=bad frame              */
    float               exc[],      /* i/o:   exc[-(PIT_MAX+L_INTERPOL)..lg] */
    float               synth[],    /* i/o:   synth[-M..lg]            */
    int                 T_out[],    /* out:   LTP lag for high band    */
    float               p_out[],    /* out:   LTP gain for high band   */
    int *               pT,         /* out:   pitch for all subframe   */
    float *             pgainT,     /* out:   pitch gain for all subfr */
    int                 pit_adj,
    float               stab_fac, /* input: stability of isf         */
    Decoder_State_Plus *st );     /* i/o :  coder memory state       */

/*---------------------------------------------------------------------*
 *              TCX routines                                           *
 *---------------------------------------------------------------------*/

void coder_tcx( float A[],       /* input: coefficients NxAz[M+1] */
    float             speech[],  /* input: speech[-M..lg] */
    float *           mem_wsp,   /* in/out: wsp memory */
    float *           mem_wsyn,  /* in/out: wsyn memory */
    float             synth[],   /* in/out: synth[-M..lg] */
    float             exc[],     /* output: exc[0..lg] */
    float             wovlp[],   /* i/o: wovlp[0..127] */
    int               ovlp_size, /* input: 0, 64 or 128 (0=acelp) */
    int               lg,        /* input: frame length */
    int               nb_bits,   /* input: number of bits allowed */
    int               prm[] );                 /* output: tcx parameters */

void decoder_tcx( int   prm[],       /* input: parameters */
    int                 nbits_AVQ[], /* input: nbits in parameters of AVQ */
    float               A[],         /* input: coefficients NxAz[M+1] */
    int                 lg,          /* input: frame length */
    int                 bad_frame[],
    float               exc[],   /* output: exc[-lg..lg] */
    float               synth[], /* in/out: synth[-M..lg] */
    Decoder_State_Plus *st );    /* i/o : coder memory state */

/*---------------------------------------------------------------------*
 *             misc                                                    *
 *---------------------------------------------------------------------*/

float segsnr(   /* return: segmential signal-to-noise ratio in dB */
    float x[],  /* input : input sequence of length n samples */
    float xe[], /* input : estimate of x */
    short n,    /* input : signal length */
    short nseg  /* input : segment length */
);

float get_gain(   /* output: codebook gain (adaptive or fixed) */
    float x[],    /* input : target signal */
    float y[],    /* input : filtered codebook excitation */
    int   L_subfr /* input : subframe size */
);

/*----------------------------------------------*
 * LPC routines.                                *
 *----------------------------------------------*/

void cos_window( float *fh, int n1, int n2 );

void q_isf_hf_new( float *isf1,  /* input : ISF in the frequency domain (0..6400) */
    float *               isf_q, /* output: quantized ISF */
    int *                 indice /* output: quantization indices (7 words) */
);
void d_isf_hf_new( int *indice, /* input: quantization indices */
    float *             isf_q,  /* output: quantized ISFs in the cosine domain */
    float *             isfold, /* input : past quantized ISF */
    int                 bfi     /* input : Bad frame indicator */
);

void q_isf_hf_16k_new( float *isf1,  /* input : ISF in the frequency domain (0..6400) */
    float *                   isf_q, /* output: quantized ISF */
    int *                     indice /* output: quantization indices (7 words) */
);
void d_isf_hf_16k_new( int *indice, /* input: quantization indices */
    float *                 isf_q,  /* output: quantized ISFs in the cosine domain */
    float *                 isfold, /* input : past quantized ISF */
    int                     bfi     /* input : Bad frame indicator */
);

int q_gain_hf_new(  /* output: indices */
    float *gain,    /* input : gain of 4 subfr */
    float *gainold, /* input : past quantized gains */
    float *gain_q   /* output: quantized gains */
);

void d_gain_hf_new( int indice,  /* input: quantization indices */
    float *             gain_q,  /* output: quantized gains */
    float *             gainold, /* input : past quantized gains */
    int                 bfi      /* input : Bad frame indicator */
);

/*---------------------------------------------------------------------*
 *                          TCX.H                                      *
 *---------------------------------------------------------------------*
 *             Prototypes of signal processing routines                *
 *---------------------------------------------------------------------*/

float AVQ_cod(      /* output: comfort noise gain factor */
    float *xri,     /* input: vector to quantize */
    int *  xriq,    /* output: quantized vector (normalized) */
    int    NB_BITS, /* input: number of bits allowed */
    int    Nsv );      /* input: number of subvector (lg=Nsv*8) */
void  AVQ_encmux( int n_pack, int *xriq, int *param, int *n_bits, int Nsv );
void  AVQ_demuxdec( int n_pack, int *param, int *n_bits, float *xriq, int Nsv, int *bfi );

int q_gain2_plus(     /* (o) : index of quantizer */
    float  code[],    /* (i) : Innovative code vector */
    int    lcode,     /* (i) : Subframe size */
    float *gain_pit,  /* (i/o): Pitch gain / Quantized pitch gain */
    float *gain_code, /* (i/o): code gain / Quantized codebook gain */
    float *coeff,     /* (i) : correlations <y1,y1>, -2<xn,y1>, */
                      /*       <y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float  mean_ener, /* (i) : mean_ener defined in open-loop (2 bits) */
    float *g0 );      /* (o)  : 'correction factor'                    */

float d_gain2_plus(      /* (o) : 'correction factor' */
    int    index,        /* (i) : index of quantizer */
    float  code[],       /* (i) : Innovative code vector */
    int    lcode,        /* (i) : Subframe size */
    float *gain_pit,     /* (o) : Quantized pitch gain */
    float *gain_code,    /* (o) : Quantized codebook gain */
    int    bfi,          /* (i) : Bad frame indicato */
    float  mean_ener,    /* (i) : mean_ener defined in open-loop (2 bits) */
    float *past_gpit,    /* (i) : past gain of pitch */
    float *past_gcode ); /* (i/o): past gain of code */

void extrapol_mod( int *bad_frame, int last_mode, int *mod );

float match_gain_6k4_new( float *AqLF, float *AqHF, int lg );

void smooth_ener_hf_new( float *HF, float *threshold, int l_frame );

void find_wsp( float A[], float speech[], /* speech[-M..lg] */
    float  wsp[],                         /* wsp[0..lg] */
    float *mem_wsp,                       /* memory */
    int    lg );

/*----------------------------------------------*
 * LPC routines.                                *
 *----------------------------------------------*/

/* new */
void init_lag_wind_hb( float bwe,    /* input : bandwidth expansion */
    float                    f_samp, /* input : sampling frequency */
    float                    wnc,    /* input : white noise correction factor */
    int                      m       /* input : order of LP filter */
);

/* new */
void lag_wind_hb( float r[], /* in/out: autocorrelations */
    int                 m    /* input : order of LP filter */
);

void isf2isp( float isf[], /* input : isf[m] normalized (range: 0<=val<=6400) */
    float           isp[], /* output: isp[m] (range: -1<=val<1) */
    int             m      /* input : LPC order */
);

/* new */
void isp2isf_hb( float isp[], /* input : isf[m] normalized (range: 0<=val<=6400) */
    float              isf[], /* output: isp[m] (range: -1<=val<1) */
    int                m      /* input : LPC order */
);

/* new */
void isp2isf_hb_16k( float isp[], /* input : isp[m] (range: -1<=val<1) */
    float                  isf[], /* output: isf[m] normalized (range: 0<=val<=5600) */
    int                    m      /* input : LPC order */
);

/* new */
void int_isp( float isp_old[], /* input : isps from past frame */
    float           isp_new[], /* input : isps from present frame */
    float           a[],       /* output: LP coefficients in both subframes */
    float           frac[],    /* input : fraction for all subframe */
    int             nb_subfr,  /* input: number of subframe */
    int             m          /* input : order of LP filter */
);

/* new */
void isf2isp_hb( float isf[], /* input : isf[m] normalized (range: 0<=val<=5600) */
    float              isp[], /* output: isp[m] (range: -1<=val<1) */
    int                m      /* input : LPC order */
);

/* new */
void isf2isp_hb_16k( float isf[], /* input : isf[m] normalized (range: 0<=val<=5600) */
    float                  isp[], /* output: isp[m] (range: -1<=val<1) */
    int                    m      /* input : LPC order */
);

void init_q_gain2( float *mem /* output :static memory (4 words) */
);

int q_gain2(          /* (o) : index of quantizer */
    float  code[],    /* (i) : Innovative code vector */
    int    lcode,     /* (i) : Subframe size */
    int    nbits,     /* (i) : number of bits (6 or 7) */
    float *gain_pit,  /* (i/o): Pitch gain / Quantized pitch gain */
    float *gain_code, /* (o) : Quantized codebook gain */
    float *coeff,     /* (i) : correlations <y1,y1>, -2<xn,y1>, */
    /* <y2,y2>, -2<xn,y2> and 2<y1,y2> */
    int    gp_clip, /* (i) : gain pitch clipping flag (1 = clipping) */
    float *mem      /* (i/o):static memory (4 words) */
);

void init_d_gain2( float *mem /* output :static memory (6 words) */
);
void d_gain2( int indice,         /* (i) : Quantization index */
    int           nbits,          /* (i) : number of bits (6 or 7) */
    float         code[],         /* (i) : Innovative code vector */
    int           lcode,          /* (i) : Subframe size */
    float *       gain_pit,       /* (o) : Quantized pitch gain */
    float *       gain_code,      /* (o) : Quantized codeebook gain */
    int           bfi,            /* (i) : Bad frame indicator */
    int           prev_bfi,       /* (i) : Previous BF indicator */
    int           state,          /* (i) : State of BFH */
    short         unusable_frame, /* (i) : UF indicator */
    short         vad_hist,       /* (i) : number of non-speech frames */
    int           i_subfr,
    float *       mem /* (i/o):static memory (6 words) */
);

void copy_coder_state(
    Coder_State_Plus *wbP, /* AMR-WB+ state struct */
    void *            st,  /* AMR-WB state struct  */
    short             sw,  /* sw=0 -> switch from WB to WB+, sw=1 -> switch from WB+ to WB */
    short             use_case_mode );

void copy_decoder_state(
    Decoder_State_Plus *wbP,
    void *              st,
    short               sw );

void  WriteHeader( EncoderConfig conf, short length, short offset, FILE *f_serial );
void  WriteBitstreamPlus( EncoderConfig conf, short length, short offset, short *serial, FILE *f_serial );
void  WriteBitstream( EncoderConfig conf, short length, short offset, unsigned char *serial, FILE *f_serial );
short ReadRawFile( short *tfi, int *bfi, DecoderConfig *conf, short *extension, short *mode, short *st_mode, short *fst, FILE *f_serial, void *serial );
int   get_nb_bits( short extension, short mode, short st_mode );
short ReadHeader( short *tfi, int *bfi, short FileFormat, short *extension, short *mode, short *st_mode, short *fst, short offset, FILE *f_serial );

#endif
