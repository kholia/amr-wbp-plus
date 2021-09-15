#ifndef int3gpp_h
#define int3gpp_h

#include "../lib_amr/typedef.h"

#define NB_QUA_GAIN7B 128
extern const Float32 E_ROM_qua_gain7b[];

void E_UTIL_residuPlus( Float32 *a, Word32 m, Float32 *x, Float32 *y, Word32 l );
void E_UTIL_synthesisPlus( Float32 a[], Word32 m, Float32 x[], Float32 y[], Word32 l, Float32 mem[], Word32 update_m );
void E_LPC_isf_reorderPlus( float *isf, float min_dist, int n );
void E_UTIL_autocorrPlus( float *x, float *r, int m, int n, float *fh );

void    E_UTIL_residu( Float32 *a, Float32 *x, Float32 *y, Word32 l );
void    E_UTIL_f_convolve( Float32 x[], Float32 h[], Float32 y[] );
void    E_UTIL_synthesis( Float32 a[], Float32 x[], Float32 y[], Word32 l, Float32 mem[], Word32 update_m );
void    E_LPC_a_isp_conversion( Float32 *a, Float32 *isp, Float32 *old_isp, Word32 m );
void    E_LPC_f_isp_a_conversion( Float32 *isp, Float32 *a, Word32 m );
void    E_LPC_lev_dur( Float32 *a, Float32 *r, Word32 m );
void    E_LPC_a_weight( Float32 *a, Float32 *ap, Float32 gamma, Word32 m );
Word16  E_UTIL_random( Word16 *seed );
void    E_GAIN_f_pitch_sharpening( Float32 *x, Word32 pit_lag );
void    E_ACELP_xh_corr( Float32 *x, Float32 *y, Float32 *h );
void    E_ACELP_xy2_corr( Float32 xn[], Float32 y1[], Float32 y2[], Float32 g_corr[] );
Float32 E_ACELP_xy1_corr( Float32 xn[], Float32 y1[], Float32 g_corr[] );
void    E_UTIL_f_preemph( Float32 *signal, Float32 mu, Word32 L, Float32 *mem );
void    E_UTIL_deemph( Float32 *signal, Float32 mu, Word32 L, Float32 *mem );
void    E_UTIL_hp50_12k8( Float32 signal[], Word32 lg, Float32 mem[] );
void    E_ACELP_codebook_target_update( Float32 *x, Float32 *x2, Float32 *y, Float32 gain );
void    E_ACELP_2t( Float32 dn[], Float32 cn[], Float32 H[], Word16 code[], Float32 y[], Word32 *index );
void    E_ACELP_4t( Float32 dn[], Float32 cn[], Float32 H[], Word16 code[], Float32 y[], Word32 nbbits, Word16 mode, Word32 _index[] );

void   D_ACELP_decode_4t( Word16 index[], Word16 nbbits, Word16 code[] );
void   E_LPC_isp_isf_conversion( Float32 isp[], Float32 isf[], Word32 m );
Word32 E_GAIN_open_loop_search( Float32 *wsp, Word32 L_min, Word32 L_max, Word32 nFrame, Word32 L_0, Float32 *gain, Float32 *hp_wsp_mem, Float32 hp_old_wsp[], UWord8 weight_flg );
Word32 E_GAIN_olag_median( Word32 prev_ol_lag, Word32 old_ol_lag[5] );
void   E_GAIN_clip_init( Float32 mem[] );
Word32 E_GAIN_clip_test( Float32 mem[] );
void   E_GAIN_clip_isf_test( Float32 isf[], Float32 mem[] );

Word32 E_GAIN_closed_loop_search( Float32 exc[], Float32 xn[], Float32 h[], Word32 t0_min, Word32 t0_max, Word32 *pit_frac, Word32 i_subfr, Word32 t0_fr2, Word32 t0_fr1 );
void   E_GAIN_adaptive_codebook_excitation( Word16 exc[], Word16 T0, Word32 frac, Word16 L_subfr );
void   E_GAIN_lp_decim2( Float32 x[], Word32 l, Float32 *mem );

#endif