#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#define L_EXTRA 96
static int short_pitch_tracker(
    float syn[], /* input:  synthesis [-PIT_MAX..L_SUBFR] */
    int   T );     /* input:  pitch period (>= PIT_MIN)     */
void init_bass_postfilter( Decoder_State_Plus *st )
{
	set_zero( st->old_synth_pf, PIT_MAX_MAX + ( 2 * L_SUBFR ) );
	set_zero( st->old_noise_pf, 2 * L_FILT );
	st->old_gain_pf[0] = 0.0;
	st->old_T_pf[0] = 64;
	st->old_gain_pf[1] = 0.0;
	st->old_T_pf[1] = 64;
	return;
}
void bass_postfilter(
    float *synth_in,  /* (i) : 12.8kHz synthesis to postfilter             */
    int *  T_sf,      /* (i) : Pitch period for all subframe (T_sf[16])    */
    float *gainT_sf,  /* (i) : Pitch gain for all subframe (gainT_sf[16])  */
    float *synth_out, /* (o) : filtered synthesis (with delay=L_SUBFR+L_FILT) */
                      /* delay = (2*L_SUBFR)+L_FILT   using EXTENSION_VA   */
    int                 pit_adj,
    Decoder_State_Plus *st ) /* i/o : decoder memory state                        */
{
	int   i, j, i_subfr, T, lg;
	float tmp, ener, gain;
	float syn_buf[PIT_MAX_MAX + ( 2 * L_SUBFR ) + L_FRAME_PLUS], *syn;
	float noise_buf[( 2 * L_FILT ) + L_SUBFR], *noise;
	int   PIT_MIN; /* Minimum pitch lag with resolution 1/4      */
	int   PIT_FR2; /* Minimum pitch lag with resolution 1/2      */
	int   PIT_FR1; /* Minimum pitch lag with resolution 1        */
	int   PIT_MAX; /* Maximum pitch lag                          */
	if( pit_adj == 0 ) {
		PIT_MIN = PIT_MIN_12k8;
		PIT_FR2 = PIT_FR2_12k8;
		PIT_FR1 = PIT_FR1_12k8;
		PIT_MAX = PIT_MAX_12k8;
	}
	else {
		i = ( ( ( pit_adj * PIT_MIN_12k8 ) + ( FSCALE_DENOM / 2 ) ) / FSCALE_DENOM ) - PIT_MIN_12k8;
		PIT_MIN = PIT_MIN_12k8 + i;
		PIT_FR2 = PIT_FR2_12k8 - i;
		PIT_FR1 = PIT_FR1_12k8;
		PIT_MAX = PIT_MAX_12k8 + ( 6 * i );
	}
	syn = syn_buf + PIT_MAX_MAX;
	noise = noise_buf + L_FILT;
	if( pit_adj == 0 ) {
		mvr2r( st->old_synth_pf, syn_buf, PIT_MAX_MAX + L_SUBFR );
		mvr2r( synth_in, syn_buf + PIT_MAX_MAX + L_SUBFR, L_FRAME_PLUS );
		mvr2r( syn_buf + L_FRAME_PLUS, st->old_synth_pf, PIT_MAX_MAX + L_SUBFR );
		mvr2r( &syn[-L_FILT], synth_out, L_FRAME_PLUS );
		T = st->old_T_pf[0];
		gain = st->old_gain_pf[0];
	}
	else {
		mvr2r( st->old_synth_pf, syn_buf, PIT_MAX_MAX + ( 2 * L_SUBFR ) );
		mvr2r( synth_in, syn_buf + PIT_MAX_MAX + ( 2 * L_SUBFR ), L_FRAME_PLUS );
		mvr2r( syn_buf + L_FRAME_PLUS, st->old_synth_pf, PIT_MAX_MAX + ( 2 * L_SUBFR ) );
		mvr2r( &syn[-L_FILT], synth_out, L_FRAME_PLUS );
	}
	for( i_subfr = 0; i_subfr < L_FRAME_PLUS; i_subfr += L_SUBFR ) {
		if( pit_adj != 0 ) {
			i = ( i_subfr / 64 ) - 2;
			if( i < 0 ) {
				T = st->old_T_pf[i + 2];
				gain = st->old_gain_pf[i + 2];
			}
			else {
				T = T_sf[i];
				gain = gainT_sf[i];
			}
		}
		if( gain > 1.0 ) {
			gain = 1.0;
		}
		if( gain < 0.0 ) {
			gain = 0.0;
		}
		T = short_pitch_tracker( &syn[i_subfr], T );
		if( pit_adj == 0 ) {
			lg = L_SUBFR + L_FRAME_PLUS - T - i_subfr;
		}
		else {
			lg = ( 2 * L_SUBFR ) + L_FRAME_PLUS - T - i_subfr;
		}
		if( lg < 0 ) {
			lg = 0;
		}
		if( lg > L_SUBFR ) {
			lg = L_SUBFR;
		}
		/* limit gain to avoid problem on burst */
		if( lg > 0 ) {
			tmp = 0.01f;
			for( i = 0; i < lg; i++ ) {
				tmp += syn[i + i_subfr] * syn[i + i_subfr];
			}
			ener = 0.01f;
			for( i = 0; i < lg; i++ ) {
				ener += syn[i + i_subfr + T] * syn[i + i_subfr + T];
			}
			tmp = (float)sqrt( tmp / ener );
			if( tmp < gain ) {
				gain = tmp;
			}
		}
		/* calculate noise based on voiced pitch */
		tmp = gain * 0.5f;
		for( i = 0; i < lg; i++ ) {
			noise_buf[i + ( 2 * L_FILT )] = tmp * ( syn[i + i_subfr] - 0.5f * syn[i + i_subfr - T] - 0.5f * syn[i + i_subfr + T] );
		}
		for( i = lg; i < L_SUBFR; i++ ) {
			noise_buf[i + ( 2 * L_FILT )] = tmp * ( syn[i + i_subfr] - syn[i + i_subfr - T] );
		}
		mvr2r( st->old_noise_pf, noise_buf, 2 * L_FILT );
		mvr2r( noise_buf + L_SUBFR, st->old_noise_pf, 2 * L_FILT );
		/* substract from voiced speech low-pass filtered noise */
		for( i = 0; i < L_SUBFR; i++ ) {
			tmp = filt_lp[0] * noise[i];
			for( j = 1; j <= L_FILT; j++ ) {
				tmp += filt_lp[j] * ( noise[i - j] + noise[i + j] );
			}
			synth_out[i + i_subfr] -= tmp;
		}
		if( pit_adj == 0 ) {
			i = i_subfr / 64;
			T = T_sf[i];
			gain = gainT_sf[i];
		}
	}
	if( pit_adj == 0 ) {
		st->old_T_pf[0] = T;
		st->old_gain_pf[0] = gain;
	}
	else {
		st->old_T_pf[0] = T_sf[NB_SUBFR - 2];
		st->old_gain_pf[0] = gainT_sf[NB_SUBFR - 2];
		st->old_T_pf[1] = T_sf[NB_SUBFR - 1];
		st->old_gain_pf[1] = gainT_sf[NB_SUBFR - 1];
	}
	return;
}
static int short_pitch_tracker(
    float syn[], /* input:  synthesis [-PIT_MAX..L_SUBFR] */
    int   T )      /* input:  pitch period (>= PIT_MIN)     */
{
	int    i, T2;
	float  tmp, corr, ener, cn;
	float *v1, *v2;
	/*----------------------------------------------------------------*
  * Test pitch/2 to avoid continuous pitch doubling                *
  * (short pitch is limited to PIT_MIN (34 = 376Hz) by the encoder *
  *----------------------------------------------------------------*/
	T2 = T >> 1;
	v1 = &syn[-L_EXTRA];
	v2 = &syn[-T2 - L_EXTRA];
	ener = 0.01f;
	for( i = 0; i < L_SUBFR + L_EXTRA; i++ ) {
		ener += v1[i] * v1[i];
	}
	corr = 0.01f;
	for( i = 0; i < L_SUBFR + L_EXTRA; i++ ) {
		corr += v1[i] * v2[i];
	}
	tmp = 0.01f;
	for( i = 0; i < L_SUBFR + L_EXTRA; i++ ) {
		tmp += v2[i] * v2[i];
	}
	/* cn = normalized correlation of pitch/2 */
	cn = corr / (float)sqrt( ener * tmp );
	if( cn > 0.95f ) {
		T = T2;
	}
	return ( T );
}
