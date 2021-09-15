#ifndef mem_h
#define mem_h

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "cnst.h"
#include "cnst_tcx_stereo.h"
#include "consts.h"
#include "nclass.h"
#include "util_stereo_x.h"
#include "wb_vad.h"

// CR3GP
typedef struct
{
	short mode;          /* AMR_WB core mode: 0..8 */
	short extension;     /* 0=AMRWB, 1=mono, 2=stereo20%, 3=stereo25% */
	short st_mode;       /* stereo mode 0..13 (not used, for ericsson cmd line?) */
	short fscale;        /* frequency scaling */
	short use_case_mode; /* use case (for AMRWB+ only) */
	short allow_dtx;     /* dtx (for AMRWB only) */
	short FileFormat;    /* File format */
	short mode_index;
	short fscale_index;
	short bc; /* Backward compatible file format*/

} EncoderConfig;

typedef struct
{
	short mode;      /* AMR_WB core mode: 0..8 */
	short extension; /* 0=AMRWB, 1=mono, 2=stereo20%, 3=stereo25% */
	short st_mode;   /* stereo mode 0..13 (not used, for ericsson cmd line?) */
	short fscale;
	long  fs;
	int   mono_dec_stereo;
	int   limiter_on;
	short FileFormat;
	short fer_sim; /* frame errasures simulation */

} DecoderConfig;

/*---------------------------------------------------------------*
 * Encoder Static RAM		  								     *
 *---------------------------------------------------------------*/

typedef struct {

	/* cod_main.c */

	float mem_decim[L_MEM_DECIM_SPLIT]; /* speech decimated filter memory */
	int   decim_frac;
	float mem_sig_in[4]; /* hp50 filter memory */
	float mem_preemph;   /* speech preemph filter memory */

	float mem_decim_hf[2 * L_FILT24k];    /* HF speech decimated filter memory */
	float old_speech_hf[L_OLD_SPEECH_ST]; /* HF old speech vector at 12.8kHz */

	/* cod_hf.c */
	float past_q_isf_hf[MHF]; /* HF past quantized isf */
	float ispold_hf[MHF];     /* HF old isp (immittance spectral pairs) */
	float ispold_q_hf[MHF];   /* HF quantized old isp */
	float old_gain;           /* HF old gain (match amplitude at 6.4kHz) */
	float mem_hf1[MHF];       /* HF memory for gain calculcation */
	float mem_hf2[MHF];       /* HF memory for gain calculcation */
	float mem_hf3[MHF];       /* HF memory for gain calculcation */

	float old_exc[PIT_MAX_MAX];

	/* Table pointers */
	const float *mean_isf_hf;
	const float *dico1_isf_hf;

} Coder_StState;

typedef struct {
	/* memory for both channels */
	Coder_StState left;
	Coder_StState right;

	/* memory for the  stereo */

	float old_chan[L_OLD_SPEECH_ST];
	float old_chan_2k[L_OLD_SPEECH_2k];
	float old_chan_hi[L_OLD_SPEECH_hi];

	float old_speech_2k[L_OLD_SPEECH_2k];
	float old_speech_hi[L_OLD_SPEECH_hi];

	float old_speech_pe[L_OLD_SPEECH_ST];
	// NMBC
	float        old_wh[HI_FILT_ORDER];
	float        old_wh_q[HI_FILT_ORDER];
	float        old_gm_gain[2];
	float        old_exc_mono[HI_FILT_ORDER];
	float        filt_energy_threshold;
	float        w_window[L_SUBFR];
	const PMSVQ *filt_hi_pmsvq;
	const PMSVQ *gain_hi_pmsvq;
	// E_STEREO_TCX
	int   mem_stereo_ovlp_size;
	float mem_stereo_ovlp[L_OVLP_2k];

	/* cod_main.c */
	NCLASSDATA *stClass;
	VadVars *   vadSt;
	short       vad_hist;

	float old_speech[L_OLD_SPEECH_ST]; /* old speech vector at 12.8kHz */
	float old_synth[M];                /* synthesis memory */

	/* cod_lf.c */
	float past_isfq[M];                       /* past isf quantizer */
	float old_wovlp[128];                     /* last tcx overlap synthesis */
	float old_d_wsp[PIT_MAX_MAX / OPL_DECIM]; /* Weighted speech vector */
	float old_exc[PIT_MAX_MAX + L_INTERPOL];  /* old excitation vector */

	float old_mem_wsyn;  /* weighted synthesis memory */
	float old_mem_w0;    /* weighted speech memory */
	float old_mem_xnq;   /* quantized target memory */
	int   old_ovlp_size; /* last tcx overlap size */

	float isfold[M];        /* old isf (frequency domain) */
	float ispold[M];        /* old isp (immittance spectral pairs) */
	float ispold_q[M];      /* quantized old isp */
	float mem_wsp;          /* wsp vector memory */
	float mem_lp_decim2[3]; /* wsp decimation filter memory */

	/* memory of open-loop LTP */

	float ada_w;
	float ol_gain;
	short ol_wght_flg;

	long int old_ol_lag[5];
	int      old_T0_med;
	float    hp_old_wsp[L_FRAME_PLUS / OPL_DECIM + ( PIT_MAX_MAX / OPL_DECIM )];
	float    hp_ol_ltp_mem[/* HP_ORDER*2 */ 3 * 2 + 1];

	//	Mem_OL mem_ol;

	/* LP analysis window */
	float window[L_WINDOW_HIGH_RATE];

	// Memory of past gain for WB+ -> WB switching
	short SwitchFlagPlusToWB;
	float mem_gain_code[4];
	short prev_mod;
} Coder_State_Plus;

/*---------------------------------------------------------------*
 * Decoder Static RAM		  								     *
 *---------------------------------------------------------------*/

typedef struct {
	/* dec_main.c */

	float mem_oversamp[L_MEM_JOIN_OVER]; /* memory for core oversampling */
	int   over_frac;
	float mem_oversamp_hf[2 * L_FILT]; /* memory for HF oversampling */

	/* dec_hf.c */
	float past_q_isf_hf[MHF];       /* HF past quantized isf */
	float past_q_isf_hf_other[MHF]; /* HF past quantized isf for the other channel when mono decoding stereo */
	float past_q_gain_hf;           /* HF past quantized gain */
	float past_q_gain_hf_other;     /* HF past quantized gain for the other channel when mono decoding stereo */
	float old_gain;                 /* HF old gain (match amplitude at 6.4kHz) */
	float ispold_hf[MHF];           /* HF old isp (immittance spectral pairs) */
	float threshold;                /* HF memory for smooth ener */
	float mem_syn_hf[MHF];          /* HF synthesis memory */

	float mem_d_tcx[D_NC + ( D_STEREO_TCX * 32 / 5 )];

	// NMBC
	float mem_d_nonc[D_NC];
	float mem_synth_hi[M];

	float mem_sig_out[4]; /* hp50 filter memory for synthesis */

	float old_synth_hf[D_BPF + L_SUBFR + L_BSP + 2 * D_NC + L_FDEL + 32 * D_STEREO_TCX / 5]; /* HF memory for synchronisation */

	float lp_amp; /* HF memory for soft exc */

	/* Table pointers */
	const float *mean_isf_hf;
	const float *dico1_isf_hf;

} Decoder_StState;

/* Memory structure for parametric stereo decoding states. */

typedef struct {
	/* memory for both channels */
	Decoder_StState left;
	Decoder_StState right;

	/* memory for parametric stereo */

	float mem_left_2k[2 * L_FDEL_2k];
	float mem_right_2k[2 * L_FDEL_2k];
	float mem_left_hi[L_FDEL];
	float mem_right_hi[L_FDEL];

	// alternative

	float my_old_synth_2k[L_FDEL_2k + D_STEREO_TCX + 2 * ( D_NC * 5 ) / 32];

	float my_old_synth_hi[2 * L_FDEL];
	float my_old_synth[2 * L_FDEL + 20];

	// NMBC
	float        old_AqLF[5 * ( M + 1 )];
	float        old_wh[HI_FILT_ORDER];
	float        old_wh2[HI_FILT_ORDER];
	float        old_exc_mono[HI_FILT_ORDER];
	float        old_gain_left[4];
	float        old_gain_right[4];
	float        old_wh_q[HI_FILT_ORDER];
	float        old_gm_gain[2];
	float        w_window[L_SUBFR];
	const PMSVQ *filt_hi_pmsvq;
	const PMSVQ *gain_hi_pmsvq;
	// E_STEREO_TCX
	int   mem_stereo_ovlp_size;
	float mem_stereo_ovlp[L_OVLP_2k];
	int   last_stereo_mode;
	float side_rms;
	float h[ECU_WIEN_ORD + 1];
	float mem_balance;
	//E_TCX_FILL
	float old_xri[L_TCX];
	/* memory for lower band (mono) */
	/* dec_main.c */
	int   last_mode;      /* last mode in previous 80ms frame */
	float mem_sig_out[4]; /* hp50 filter memory for synthesis */
	float mem_deemph;     /* speech deemph filter memory      */
	/* dec_lf.c */
	int   prev_lpc_lost; /* previous lpc is lost when = 1 */
	float old_synth[M];  /* synthesis memory */

	float old_exc[PIT_MAX_MAX + L_INTERPOL]; /* old excitation vector (>20ms) */

	float isfold[M];                      /* old isf (frequency domain) */
	float ispold[M];                      /* old isp (immittance spectral pairs) */
	float past_isfq[M];                   /* past isf quantizer */
	float wovlp[128];                     /* last weighted synthesis for overlap */
	int   ovlp_size;                      /* overlap size */
	float isf_buf[L_MEANBUF * ( M + 1 )]; /* old isf (for frame recovery) */
	int   old_T0;                         /* old pitch value (for frame recovery) */
	int   old_T0_frac;                    /* old pitch value (for frame recovery) */
	short seed_ace;                       /* seed memory (for random function) */
	float mem_wsyn;                       /* TCX synthesis memory */
	short seed_tcx;                       /* seed memory (for random function) */
	float wsyn_rms;                       /* rms value of weighted synthesis */
	float past_gpit;                      /* past gain of pitch (for frame recovery) */
	float past_gcode;                     /* past gain of code (for frame recovery) */
	int   pitch_tcx;                      /* for bfi on tcx20 */
	float gc_threshold;

	/* bass_pf.c */
	float old_synth_pf[PIT_MAX_MAX + ( 2 * L_SUBFR )]; /* bass post-filter: old synthesis  */
	float old_noise_pf[2 * L_FILT];                    /* bass post-filter: noise memory   */
	int   old_T_pf[2];                                 /* bass post-filter: old pitch      */
	float old_gain_pf[2];                              /* bass post-filter: old pitch gain */
	/* Table pointers */
	float *mean_isf_hf;
	float *dico1_isf_hf;

	// For WB <-> WB+ switching
	float mem_gain_code[4];
	float mem_lpc_hf[MHF + 1];
	float mem_gain_hf;
	short ramp_state;

} Decoder_State_Plus;

#endif /* mem_h */
