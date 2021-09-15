/* Stereo utilities */

#ifndef STEREO_UTIL_H
#define STEREO_UTIL_H
#include "cnst.h"
#include "mem.h" // needed for definition of "Coder_State_Plus"
#include "options.h"

float my_max( float x, float y );
float my_min( float x, float y );

void crosscorrVA( float *x, /* input : input signal */
    float *              y,
    float *              r, /* output: autocorrelations vector */
    int                  m, /* input : order of LP filter */
    int                  n, /* input : window size */
    float *              fh /* input : analysis window */
);

void fir_filt( float *h, /* input : filter coefficients */
    int               m, /* input : order of filter */
    float *           x, /* input : input signal (usually mono) */
    float *           y, /* output: output signal (usually stereo) */
    int               l  /* input : size of filtering */
);

void int_fir( float fir_old[], /* input : FIR from past frame */
    float           fir_new[], /* input : FIR from present frame */
    float           h[],       /* output: FIR coefficients in both subframes */
    int             nb_subfr,  /* input: number of subframe */
    int             m          /* input : order of FIR filter */
);

void mix_ch( float *ch_left,  /* input: samples from left channel */
    float *         ch_right, /* input: samples from right channel */
    float *         ch_sum,   /* output: mixed mono signal */
    int             n,        /* input: length of frame */
    float           gA,       /* input: weight factor for left channel */
    float           gB        /* input: weight factor for right channel */
);

#ifndef CLEAN_UP_VCC_WARN
#ifdef COSINE_INTERP
void init_csbi_window( void );
#endif
#ifdef LAGW
void init_lag_wind_su( float lag_window[], float bwe, /* input : bandwidth expansion */
    float f_samp,                                     /* input : sampling frequency */
    float wnc,                                        /* input : white noise correction factor */
    int   m                                           /* input : order of LP filter */
);
#endif

#ifdef LAGW
void lag_wind_s( float lag_window[], float rcc[], float rca[], float rcb[], int m );
void lag_wind_s2( float lag_window[], float rcc[], float rca[], int m );

#else
void init_lag_wind_s( float bwe, float f_samp, float wnc, int m );

void lag_wind_s( float r[], int m );

#endif
#endif

int chol( int order, /* (i) : Order */
    float *   a,
    float *   c,
    float *   b );

/*-------------------------------------------------------------------*
 * Function  band_split												 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~										 *
 *   ->	Routine for splitting one frame into one high-frequency part *
 *		and one low-frequency part.									 *
 *-------------------------------------------------------------------*/

void band_split( float *x_fb,  /* (i) : Full-band samples */
    float *             x_hb,  /* (o) : Hi-band samples */
    float *             x_lb,  /* (o) : Low-band samples */
    float *             h,     /* (i) : Filter coefficients */
    int                 ncoeff /* (i) : Number of filter coefficients */
);

void band_split_2k( float sig[], float sig_2k[], float sig_hi[], int lg, float mem_sig[], float mem_sig_2k[] );

void band_split_taligned_2k( float sig[], float sig_2k[], float sig_hi[], int lg );

/*-------------------------------------------------------------------*
 * Function  band_join												 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~										 *
 *   ->	Routine for joining one high-frequency part					 *
 *		and one low-frequency part into one frame.					 *
 *-------------------------------------------------------------------*/

void band_join_2k( float sig[], float sig_2k[], float sig_hi[], int lg );

void gen_chan( float *ac,     /* (i) : filter coefficients */
    float *           c,      /* (i) : Mixed channel */
    float *           a,      /* (o) : Left channel */
    float *           b,      /* (o) : Right channel */
    int               band_2k /* (i) : Operation in low band */
);

void write_wave_header( FILE *fp, int srate, int nCh, int size );

/* Mix channels with prefiltering for sum channel*/
void mix_ch_pf( float *ch_left,   /* input: samples from left channel */
    float *            ch_right,  /* input: samples from right channel */
    float *            ch_sum,    /* output: mixed mono signal */
    float *            ch_sum_pf, /* output: mixed mono signal, prefiltered */
    Coder_State_Plus * st );       /* i/o : coder memory state */

#endif
