#ifndef table_decl_h
#define table_decl_h

#include "consts.h"

/* RE8 lattice quantiser tables */
extern const int          tab_pow2[8];
extern const int          tab_factorial[8];
extern const int          Ia[NB_LEADER];
extern const int          Ds[NB_LDSIGN];
extern const unsigned int Is[NB_LDSIGN];
extern const int          Da[][8], Ns[], A3[], A4[];
extern const unsigned int I3[], I4[];
extern const int          Da_nq[];
extern const int          Da_pos[], Da_nb[], Da_id[];

/* LF ISF quantiser codebooks */
extern const float isf_init[M];
extern const float mean_isf[16];
extern const float dico1_isf[SIZE_BK1 * 9];
extern const float dico2_isf[SIZE_BK2 * 7];
extern const float dico21_isf[SIZE_BK21 * 3];
extern const float dico22_isf[SIZE_BK22 * 3];
extern const float dico23_isf[SIZE_BK23 * 3];
extern const float dico24_isf[SIZE_BK24 * 3];
extern const float dico25_isf[SIZE_BK25 * 4];
extern const float dico21_isf_36b[SIZE_BK21_36b * 5];
extern const float dico22_isf_36b[SIZE_BK22_36b * 4];
extern const float dico23_isf_36b[SIZE_BK23_36b * 7];

/* HF ISF quantiser codebooks */
extern const float isf_init_hf[M];
extern const float mean_isf_hf_low_rate[Q_ISF_ORDER];
extern const float dico1_isf_hf_low_rate[SIZE_BK1_HF * Q_ISF_ORDER];
extern const float dico2_isf_hf[SIZE_BK2_HF * Q_ISF_ORDER];
extern const float dico1_isf_hf_low_rate[SIZE_BK1_HF * Q_ISF_ORDER];

extern const float mean_isf_hf_12k8[Q_ISF_ORDER];
extern const float dico1_isf_hf_12k8[SIZE_BK1_HF * Q_ISF_ORDER];

/* HF gain quantiser codebook */
extern const float dico_gain_hf[SIZE_BK_HF * Q_GN_ORDER];

/* Resampling filters */
extern const float filter_8k[61];
extern const float filter_28k8[109];
extern const float filter_28k8_hf[109];
extern const float filter_32k[61];
extern const float filter_32k_hf[61];
extern const float filter_32k_7k[61];
extern const float filter_48k[185];
extern const float filter_48k_hf[185];
extern const float filter_48k_hf_high_rate[185];

/* FFT sin and cos tables */
extern const float t_sin[N_MAX];
extern const float t_cos[N_MAX];

/* For phase random generator */
extern const float sin20[20];

/* For pitch predictor */
extern const float inter4_2[PIT_FIR_SIZE2];

/* For bass postfilter */
extern const float filt_lp[1 + L_FILT];

/* 2k LP filter sed by stereo code */
extern const float filter_2k[N_COEFF_F2K];

extern const float lag_window[17];
extern const float sin20[20];

/* for open-loop classifier */
extern const float VadFiltBandFreqs[COMPLEN2];
extern const float bw[COMPLEN];
extern const float lwg[LPHAVELEN];

/* stereo */
extern const int nprm_stereo_hi_x[4];

/* Bit allocation tables for stereo */
extern const int StereoNbits[16];

extern const float filter_2k[N_COEFF_F2K];

extern const float filt_lp[1 + L_FILT];

/* High stereo codebooks */
// average filter
extern const float cb_filt_hi_mean[HI_FILT_ORDER];
extern const float filt_hi_mscb_4a[SIZE_FILT_HI_MSVQ_4A][HI_FILT_ORDER];
extern const float filt_hi_mscb_7a[SIZE_FILT_HI_MSVQ_7A][HI_FILT_ORDER];
extern const float filt_hi_mscb_7b[SIZE_FILT_HI_MSVQ_7B][HI_FILT_ORDER];

// average gain vector
extern const float cb_gain_hi_mean[HI_GAIN_ORDER];
extern const float gain_hi_mscb_2a[SIZE_GAIN_HI_MSVQ_2A][HI_GAIN_ORDER];
extern const float gain_hi_mscb_5a[SIZE_GAIN_HI_MSVQ_5A][HI_GAIN_ORDER];

// code books
// 4 bit VQ (filter)
extern const int    size_filt_hi_msvq_4[NSTAGES_FILT_HI_MSVQ4];
extern const float *cbs_filt_hi_msvq_4[NSTAGES_FILT_HI_MSVQ4];
extern const PMSVQ  filt_hi_pmsvq4;
// 4+3 bit msvq (filter)
extern const int    size_filt_hi_msvq_7[NSTAGES_FILT_HI_MSVQ7];
extern const float *cbs_filt_hi_msvq_7[NSTAGES_FILT_HI_MSVQ7];
extern const PMSVQ  filt_hi_pmsvq7;
// 2 bit VQ (gain)
extern const int    size_gain_hi_msvq_2[NSTAGES_GAIN_HI_MSVQ2];
extern const float *cbs_gain_hi_msvq_2[NSTAGES_GAIN_HI_MSVQ2];
extern const PMSVQ  gain_hi_pmsvq2;
// 5 bit VQ (gain)
extern const int    size_gain_hi_msvq_5[NSTAGES_GAIN_HI_MSVQ5];
extern const float *cbs_gain_hi_msvq_5[NSTAGES_GAIN_HI_MSVQ5];
extern const PMSVQ  gain_hi_pmsvq5;

extern const float gain_hf_ramp[64];

extern const float inter2_coef[12];

extern const float filter_LP12[1 + ( L_FILT_OVER_FS * 12 )];
#ifdef FILTER_44kHz
extern const float filter_LP165[1 + ( ( 12 + 1 ) * 165 )];
#endif
#ifdef FILTER_48kHz
extern const float filter_LP180[1 + ( ( 12 + 1 ) * 180 )];
#endif
// conveniant tables, mode index and isf index as defined in TS26.290

extern const short miMode[];
extern const short isfIndex[];
extern const short StereoRate[];
extern const short MonoRate[];

#endif
