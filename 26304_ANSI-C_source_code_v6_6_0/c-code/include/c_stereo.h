
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

#ifdef E_STEREO_2
void coder_stereo( float speech_hi[], /* (i) : Mixed channel, hi */
    float                chan_hi[],   /* (i) : Right channel, hi */
    float                speech_2k[], /* (i) : Mixed channel, lo */
    float                chan_2k[],   /* (i) : Right channel, lo */
    int                  param[],     /* (o) : Encoded parameters */
    int                  brMode,
    short                d_tcx_serial[],
    Coder_State_Plus *   st ); /* (i/o): Encoder states */
#else
void coder_stereo( float side_speech[], /* (i) : L-R channel part (fs=12.8kHz) */
    float                speech[],      /* (i) : right channel part (fs=12.8kHz) */
    float                synth[],
    short                param[], /* (o) : parameters */
    float                window[],
    Coder_State_Plus *   st ); /* i/o : coder memory state */
#endif

/*-----------------------------------------------------------------*
 * Function  init_coder_stereo                                      *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~                                      *
 *   ->Initialization of variables for the stereo coder.           *
 *-----------------------------------------------------------------*/
void init_coder_stereo( Coder_State_Plus *st ); /* (i/o) : Encoder states. */

/*-----------------------------------------------------------------*
 * Function  end_prm_stereo									       *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~                                      *
 *   ->Encoding of stereo parameters.					           *
 *-----------------------------------------------------------------*/

void enc_prm_stereo( int param[],   /* (i) : parameters */
    short                serial[],  /* (o) : serial bits stream */
    int                  nbits_pack /* (i) : number of bits per packet of 20ms */
);

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

/*-------------------------------------------------------------------*
 * Function  lo_filt_design											 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~										 *
 *   ->	Design of Wiener filters for left/right channel for low-	 *
 *		frequency band (below 1 kHz). 								 *
 *-------------------------------------------------------------------*/

void lo_filt_design( float *spR,     /* (i) : Speech input (R) */
    float *                 spM,     /* (i) : Speech input (M) */
    float *                 filt_lo, /* (o) : Filter coeff output */
    Coder_State_Plus *      st       /* (i/o): Filter states */
);

/*-------------------------------------------------------------------*
 * Function  hi_filt_design											 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~										 *
 *   ->	Design of Wiener filters for left/right channel for hi-	 *
 *		frequency band (above 1 kHz). 								 *
 *-------------------------------------------------------------------*/

void hi_filt_design( float *spR,     /* (i) : Speech input (R) */
    float *                 spM,     /* (i) : Speech input (M) */
    float *                 filt_hi, /* (o) : Filter coeff output */
    Coder_State_Plus *      st       /* (i/o): Filter states */
);
