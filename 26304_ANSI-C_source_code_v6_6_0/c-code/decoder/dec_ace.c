#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
    int                 pit_adj,    /* input: pitch adjust */
    float               stab_fac,   /* input: stability of isf         */
    Decoder_State_Plus *st )        /* i/o :  coder memory state       */
{
	int    i, i_subfr, select;
	int    T0, T0_frac, index, pit_flag, T0_min, T0_max;
	float  tmp, gain_pit, gain_code, voice_fac, ener, mean_ener_code, ener_wsyn;
	float  code[L_SUBFR], buf[M + L_OVLP];
	float  fac, exc2[L_SUBFR];
	float *p_A, Ap[1 + M];
	short  code3GPP[L_SUBFR], prm3GPP[10];
	int    PIT_MIN; /* Minimum pitch lag with resolution 1/4      */
	int    PIT_FR2; /* Minimum pitch lag with resolution 1/2      */
	int    PIT_FR1; /* Minimum pitch lag with resolution 1        */
	int    PIT_MAX; /* Maximum pitch lag                          */
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
	/*------------------------------------------------------------------------*
  * - decode mean_ener_code for gain decoder (d_gain2.c)                   *
  *------------------------------------------------------------------------*/
	index = *prm++;
	/* decode mean energy with 2 bits : 18, 30, 42 or 54 dB */
	if( !bfi ) {
		mean_ener_code = ( ( (float)index ) * 12.0f ) + 18.0f;
	}
	/*------------------------------------------------------------------------*
  *          Loop for every subframe in the analysis frame                 *
  *------------------------------------------------------------------------*
  *  To find the pitch and innovation parameters. The subframe size is     *
  *  L_SUBFR and the loop is repeated L_ACELP/L_SUBFR times.               *
  *     - compute impulse response of weighted synthesis filter (h1[])     *
  *     - compute the target signal for pitch search                       *
  *     - find the closed-loop pitch parameters                            *
  *     - encode the pitch dealy                                           *
  *     - update the impulse response h1[] by including fixed-gain pitch   *
  *     - find target vector for codebook search                           *
  *     - correlation between target vector and impulse response           *
  *     - codebook search                                                  *
  *     - encode codebook address                                          *
  *     - VQ of pitch and codebook gains                                   *
  *     - find synthesis speech                                            *
  *     - update states of weighting filter                                *
  *------------------------------------------------------------------------*/
	ener_wsyn = 0.01f;
	p_A = A;
	for( i_subfr = 0; i_subfr < lg; i_subfr += L_SUBFR ) {
		pit_flag = i_subfr;
		if( i_subfr == ( 2 * L_SUBFR ) ) {
			pit_flag = 0;
		}
		index = *prm++;
		/*-------------------------------------------------*
     * - Decode pitch lag                              *
     *-------------------------------------------------*/
		if( bfi ) /* if frame erasure */
		{
			/* Lag indices received also in case of BFI, so that
         the parameter pointer stays in sync. */
			st->old_T0_frac += 1; /* use last delay incremented by 1/4 */
			if( st->old_T0_frac > 3 ) {
				st->old_T0_frac -= 4;
				( st->old_T0 )++;
			}
			if( st->old_T0 >= PIT_MAX ) {
				st->old_T0 = PIT_MAX - 5;
			}
			T0 = st->old_T0;
			T0_frac = st->old_T0_frac;
		}
		else {
			if( pit_flag == 0 ) {
				if( index < ( PIT_FR2 - PIT_MIN ) * 4 ) {
					T0 = PIT_MIN + ( index / 4 );
					T0_frac = index - ( T0 - PIT_MIN ) * 4;
				}
				else if( index < ( ( PIT_FR2 - PIT_MIN ) * 4 + ( PIT_FR1 - PIT_FR2 ) * 2 ) ) {
					index -= ( PIT_FR2 - PIT_MIN ) * 4;
					T0 = PIT_FR2 + ( index / 2 );
					T0_frac = index - ( T0 - PIT_FR2 ) * 2;
					T0_frac *= 2;
				}
				else {
					T0 = index + PIT_FR1 - ( ( PIT_FR2 - PIT_MIN ) * 4 ) - ( ( PIT_FR1 - PIT_FR2 ) * 2 );
					T0_frac = 0;
				}
				/* find T0_min and T0_max for subframe 2 or 4 */
				T0_min = T0 - 8;
				if( T0_min < PIT_MIN ) {
					T0_min = PIT_MIN;
				}
				T0_max = T0_min + 15;
				if( T0_max > PIT_MAX ) {
					T0_max = PIT_MAX;
					T0_min = T0_max - 15;
				}
			}
			else /* if subframe 2 or 4 */
			{
				T0 = T0_min + index / 4;
				T0_frac = index - ( T0 - T0_min ) * 4;
			}
		}
		/*-------------------------------------------------*
     * - Find the pitch gain, the interpolation filter *
     *   and the adaptive codebook vector.             *
     *-------------------------------------------------*/
		pred_lt4( &exc[i_subfr], T0, T0_frac, L_SUBFR + 1 );
		select = *prm++;
		if( bfi ) {
			select = 1;
		}
		if( select == 0 ) {
			/* find pitch excitation with lp filter */
			for( i = 0; i < L_SUBFR; i++ ) {
				code[i] = (float)( 0.18 * exc[i - 1 + i_subfr] + 0.64 * exc[i + i_subfr] + 0.18 * exc[i + 1 + i_subfr] );
			}
			mvr2r( code, &exc[i_subfr], L_SUBFR );
		}
		/*-------------------------------------------------------*
    * - Decode innovative codebook.                         *
    * - Add the fixed-gain pitch contribution to code[].    *
    *-------------------------------------------------------*/
		{
			int g;
			for( g = 0; g < 8; g++ ) {
				prm3GPP[g] = (short)prm[g];
			}
		}
		if( codec_mode == MODE_9k6 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 20, code3GPP );
			}
			prm += 4;
		}
		else if( codec_mode == MODE_11k2 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 28, code3GPP );
			}
			prm += 4;
		}
		else if( codec_mode == MODE_12k8 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 36, code3GPP );
			}
			prm += 4;
		}
		else if( codec_mode == MODE_14k4 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 44, code3GPP );
			}
			prm += 4;
		}
		else if( codec_mode == MODE_16k ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 52, code3GPP );
			}
			prm += 4;
		}
		else if( codec_mode == MODE_18k4 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 64, code3GPP );
			}
			prm += 8;
		}
		else if( codec_mode == MODE_20k ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 72, code3GPP );
			}
			prm += 8;
		}
		else if( codec_mode == MODE_23k2 ) {
			if( !bfi ) {
				D_ACELP_decode_4t( prm3GPP, 88, code3GPP );
			}
			prm += 8;
		}
		else {
			printf( "invalid mode for acelp frame!\n" );
			exit( 0 );
		}
		{
			int g;
			for( g = 0; g < L_SUBFR; g++ ) {
				code[g] = (float)( code3GPP[g] / 512 );
			}
		}
		if( bfi ) {
			/* the innovative code doesn't need to be scaled (see D_gain2) */
			for( i = 0; i < L_SUBFR; i++ ) {
				code[i] = (float)E_UTIL_random( &( st->seed_ace ) );
			}
		}
		/*-------------------------------------------------------*
    * - Add the fixed-gain pitch contribution to code[].    *
    *-------------------------------------------------------*/
		tmp = 0.0;
		E_UTIL_f_preemph( code, TILT_CODE, L_SUBFR, &tmp );
		i = T0;
		if( T0_frac > 2 ) {
			i++;
		}
		E_GAIN_f_pitch_sharpening( code, i );
		/*-------------------------------------------------*
    * - Decode codebooks gains.                       *
    *-------------------------------------------------*/
		index = *prm++;
		st->mem_gain_code[i_subfr / L_SUBFR] = d_gain2_plus( index, code, L_SUBFR, &gain_pit, &gain_code, bfi, mean_ener_code, &( st->past_gpit ), &( st->past_gcode ) );
		/*----------------------------------------------------------*
    * Update parameters for the next subframe.                 *
    * - tilt of code: 0.0 (unvoiced) to 0.5 (voiced)           *
    *----------------------------------------------------------*/
		/* energy of pitch excitation */
		ener = 0.0;
		for( i = 0; i < L_SUBFR; i++ ) {
			ener += exc[i + i_subfr] * exc[i + i_subfr];
		}
		ener *= ( gain_pit * gain_pit );
		/* energy of innovative code excitation */
		tmp = 0.0;
		for( i = 0; i < L_SUBFR; i++ ) {
			tmp += code[i] * code[i];
		}
		tmp *= gain_code * gain_code;
		/* find voice factor (1=voiced, -1=unvoiced) */
		voice_fac = (float)( ( ener - tmp ) / ( ener + tmp + 0.01 ) );
		/*-------------------------------------------------------*
    * - Find the total excitation.                          *
    *-------------------------------------------------------*/
		for( i = 0; i < L_SUBFR; i++ ) {
			exc2[i] = gain_pit * exc[i + i_subfr];
		}
		for( i = 0; i < L_SUBFR; i++ ) {
			exc[i + i_subfr] = gain_pit * exc[i + i_subfr] + gain_code * code[i];
		}
		/*-------------------------------------------------------*
    * - Output pitch parameters for bass post-filter        *
    *-------------------------------------------------------*/
		i = T0;
		if( T0_frac > 2 ) {
			i++;
		}
		if( i > PIT_MAX ) {
			i = PIT_MAX;
		}
		*pT++ = i;
		*pgainT++ = gain_pit;
		if( pit_adj == 0 ) {
			E_UTIL_synthesis( p_A, &exc[i_subfr], &synth[i_subfr], L_SUBFR, &synth[i_subfr - M], 0 );
		}
		else {
			/*------------------------------------------------------------*
      * noise enhancer                                             *
      * ~~~~~~~~~~~~~~                                             *
      * - Enhance excitation on noise. (modify gain of code)       *
      *   If signal is noisy and LPC filter is stable, move gain   *
      *   of code 1.5 dB toward gain of code threshold.            *
      *   This decrease by 3 dB noise energy variation.            *
      *------------------------------------------------------------*/
			tmp = (float)( 0.5 * ( 1.0 - voice_fac ) ); /* 1=unvoiced, 0=voiced */
			fac = stab_fac * tmp;
			tmp = gain_code;
			if( tmp < st->gc_threshold ) {
				tmp = (float)( tmp * 1.19 );
				if( tmp > st->gc_threshold ) {
					tmp = st->gc_threshold;
				}
			}
			else {
				tmp = (float)( tmp / 1.19 );
				if( tmp < st->gc_threshold ) {
					tmp = st->gc_threshold;
				}
			}
			st->gc_threshold = tmp;
			gain_code = (float)( ( fac * tmp ) + ( ( 1.0 - fac ) * gain_code ) );
			for( i = 0; i < L_SUBFR; i++ ) {
				code[i] *= gain_code;
			}
			/*------------------------------------------------------------*
      * pitch enhancer                                             *
      * ~~~~~~~~~~~~~~                                             *
      * - Enhance excitation on voice. (HP filtering of code)      *
      *   On voiced signal, filtering of code by a smooth fir HP   *
      *   filter to decrease energy of code in low frequency.      *
      *------------------------------------------------------------*/
			tmp = (float)( 0.125 * ( 1.0 + voice_fac ) ); /* 0.25=voiced, 0=unvoiced */
			exc2[0] += code[0] - ( tmp * code[1] );
			for( i = 1; i < L_SUBFR - 1; i++ ) {
				exc2[i] += code[i] - ( tmp * code[i - 1] ) - ( tmp * code[i + 1] );
			}
			exc2[L_SUBFR - 1] += code[L_SUBFR - 1] - ( tmp * code[L_SUBFR - 2] );
			E_UTIL_synthesis( p_A, exc2, &synth[i_subfr], L_SUBFR, &synth[i_subfr - M], 0 );
		}
		/* find weighted synthesis energy (for frame recovery on TCX mode) */
		E_LPC_a_weight( p_A, Ap, GAMMA1, M );
		E_UTIL_residu( Ap, &synth[i_subfr], code, L_SUBFR );
		E_UTIL_deemph( code, TILT_FAC, L_SUBFR, &( st->mem_wsyn ) );
		for( i = 0; i < L_SUBFR; i++ ) {
			ener_wsyn += code[i] * code[i];
		}
		p_A += ( M + 1 );
	} /* end of subframe loop */
	/* RMS of weighted synthesis (for frame recovery on TCX mode) */
	st->wsyn_rms = (float)sqrt( ener_wsyn / (float)lg );
	/*----------------------------------------------------------*
  * find 10ms ZIR in weighted domain for next tcx frame      *
  *----------------------------------------------------------*/
	mvr2r( &synth[lg - M], buf, M );
	set_zero( buf + M, L_OVLP );
	E_UTIL_synthesis( p_A, buf + M, buf + M, L_OVLP, buf, 0 );
	E_LPC_a_weight( p_A, Ap, GAMMA1, M ); /* wAi of tcx is quantized */
	E_UTIL_residu( Ap, buf + M, st->wovlp, L_OVLP );
	tmp = st->mem_wsyn;
	E_UTIL_deemph( st->wovlp, TILT_FAC, L_OVLP, &tmp );
	for( i = 1; i < ( L_OVLP / 2 ); i++ ) {
		st->wovlp[L_OVLP - i] *= ( (float)i ) / ( (float)( L_OVLP / 2 ) );
	}
	st->ovlp_size = 0; /* indicate ACELP frame to TCX */
	/* update pitch value for bfi procedure */
	st->old_T0_frac = T0_frac;
	st->old_T0 = T0;
	return;
}
