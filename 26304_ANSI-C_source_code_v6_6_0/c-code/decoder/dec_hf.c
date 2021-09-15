#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
 * Funtion  init_decoder_hf                                        *
 * ~~~~~~~~~~~~~~~~~~~~~~~~                                        *
 *                                                                 *
 *   ->Initialization of variables for the decoder section.        *
 *     (band 6400Hz to 10800Hz).                                   *
 *                                                                 *
 *-----------------------------------------------------------------*/
void init_decoder_hf( Decoder_StState *st )
{
	int i;
	/* Static vectors to zero */
	set_zero( st->mem_syn_hf, MHF );
	set_zero( st->past_q_isf_hf, MHF );
	set_zero( st->past_q_isf_hf_other, MHF );
	st->past_q_gain_hf = 0.0;
	st->past_q_gain_hf_other = 0.0;
	st->old_gain = 1.0;
	st->threshold = 0.0;
	/* Initialize the ISPs */
	for( i = 0; i < MHF - 1; i++ ) {
		st->ispold_hf[i] = (float)cos( 3.141592654 * (float)( i + 1 ) / (float)MHF );
	}
	st->ispold_hf[MHF - 1] = 0.045f;
	st->lp_amp = 0.0;
	return;
}
/*-----------------------------------------------------------------*
 * Funtion decoder_lf                                              *
 * ~~~~~~~~~~~~~~~~~~                                              *
 *   ->Principle decoder routine (band 6400Hz to 10800Hz).         *
 *                                                                 *
 *-----------------------------------------------------------------*/
void decoder_hf(
    int              mod[],           /* (i)  : mode for each 20ms frame (mode[4]     */
    int              param[],         /* (i)  : parameters                            */
    int              param_other[],   /* (i)  : parameters for the right channel in case of mono decoding with stereo bitstream */
    int              mono_dec_stereo, /* 1=Mono decoding with stereo bitstream */
    int              bad_frame[],     /* (i)  : for each frame (bad_frame[4])         */
    float            AqLF[],          /* (i)  : decoded coefficients (AdLF[16])       */
    float            exc[],           /* (i)  : decoded excitation                    */
    float            synth_hf[],      /* (o)  : decoded synthesis                     */
    float            mem_lpc_hf[],
    float *          mem_gain_hf,
    short *          ramp_state,
    Decoder_StState *st ) /* i/o : decoder memory state                */
{
	/* LPC coefficients */
	float  ispnew[MHF];
	float  isfnew[MHF];
	float  Aq[( NB_SUBFR + 1 ) * ( MHF + 1 )]; /* A(z) for all subframes     */
	float *p_Aq;
	float  gain4[4], gain_hf[NB_SUBFR], gain_hf_other[NB_SUBFR];
	float  HF[L_SUBFR];
	/* Scalars */
	int   i, k, mode, bfi, sf, nsf, *prm, *prm_other;
	int   i_subfr, index;
	float gain, tmp;
	int * prm_in_other;
	float gdiff;
	float thr_lo, thr_hi;
	float gain_sum, gain_sum_other;
	float gain_tmp[NB_SUBFR];
	float isfnew_other[MHF];
	/*----------- DECODE AND SYNTHESIZE UPPER-BAND SIGNAL (synth_hf[])
                AND ADD IT TO UPSAMPLED LOWER-BAND SIGNAL (synthesis) ----------
    bitstream pointer = ptr
  coding modes (coder type=ACELP,TCX and frame length=20,40,80 ms) for upper-band are copied from lower-band
  HF decoding depends on bit allocation to HFs (nbits_hf) :
     nbits_hf=0:
       the lower-band signal post-filtered and upsampled is the actual synthesis
       (no HF decoding and synthesis)
     nbits_hf=4*16:
       the lower-band excitation is randomized, folded, and shaped first in time-domain
       with subframe gains then in frequency domain with a LP filter (plus energy smoothing)
   */
	/* compute Ai (Aq) and correction scale factor (gain) for each subframe */
	for( k = 0; k < NB_DIV; k++ ) {
		bfi = bad_frame[k];
		mode = mod[k];
		if( mode == 0 ) {
			nsf = 4;
		}
		else {
			nsf = 2 << mode;
		}
		/* set pointer to parameters */
		prm = param + ( k * NPRM_BWE_DIV );
		if( mono_dec_stereo == 1 ) {
			prm_other = param_other + ( k * NPRM_BWE_DIV );
		}
		if( mono_dec_stereo == 1 ) {
			prm_in_other = prm_other;
			gain_sum = 0.0f;
			gain_sum_other = 0.0f;
		}
		/* decode ISFs and convert ISFs to cosine domain */
		d_isf_hf( prm, isfnew, st->past_q_isf_hf, bfi, st->mean_isf_hf, st->dico1_isf_hf );
		if( mono_dec_stereo == 0 ) {
			isf2isp( isfnew, ispnew, MHF );
			/* interpolate Ai in ISP domain (Aq) */
			int_lpc_np1( st->ispold_hf, ispnew, Aq, nsf, MHF );
			mvr2r( ispnew, st->ispold_hf, MHF );
			/* compute gain correction factor to match amplitude at 6.4kHz */
			gain = match_gain_6k4( &AqLF[( ( k * 4 ) + nsf ) * ( M + 1 )], &Aq[nsf * ( MHF + 1 )] );
			/* interpolate per subframe */
			int_gain( st->old_gain, gain, gain_hf, nsf );
			st->old_gain = gain;
		}
		/* decode 4 gains */
		d_gain_hf( prm[2], gain4, &( st->past_q_gain_hf ), bfi );
		prm += 3;
		/* decode gain correction offsets (in dB) for each subframe */
		if( mode == 3 ) { /* 80+10-ms TCX */
			/* decode gain offset: -10.5, -7.5, -4.5, -1.5, +1.5, +4.5, +7.5, +10.5 dB */
			for( i = 0; i < 16; i++ ) {
				index = *prm++;
				tmp = ( 3.0f * ( (float)index ) ) - 10.5f;
				if( bfi ) {
					tmp = 0.0;
				}
				if( bad_frame[k + 1] && ( i < 8 ) ) {
					tmp = 0.0;
				}
				if( bad_frame[k + 2] && ( i >= 8 ) ) {
					tmp = 0.0;
				}
				if( mono_dec_stereo == 1 ) {
					gain_hf[i] = tmp + gain4[( i >> 2 )];
					gain_sum += gain_hf[i];
				}
				else {
					gain_hf[i] += tmp + gain4[( i >> 2 )];
				}
			}
		}
		else if( mode == 2 ) { /* 40+5-ms TCX */
			/* decode gain offset: -4.5, -1.5, +1.5, +4.5 dB */
			for( i = 0; i < 8; i++ ) {
				index = *prm++;
				tmp = ( 3.0f * ( (float)index ) ) - 4.5f;
				if( bfi || bad_frame[k + 1] ) {
					tmp = 0.0;
				}
				if( mono_dec_stereo == 1 ) {
					gain_hf[i] = tmp + gain4[( i >> 1 )];
					gain_sum += gain_hf[i];
				}
				else {
					gain_hf[i] += tmp + gain4[( i >> 1 )];
				}
			}
		}
		else { /* ACELP or 20+2.5-ms TCX */
			if( mono_dec_stereo == 1 ) {
				for( i = 0; i < 4; i++ ) {
					gain_hf[i] = gain4[i];
					gain_sum += gain_hf[i];
				}
			}
			else {
				for( i = 0; i < 4; i++ ) {
					gain_hf[i] += gain4[i];
				}
			}
		}
		if( mono_dec_stereo == 1 ) {
			/* decode 4 gains */
			d_gain_hf( prm_other[2], gain4, &( st->past_q_gain_hf_other ), bfi );
			prm_other += 3;
			/* decode gain correction offsets (in dB) for each subframe */
			if( mode == 3 ) { /* 80+10-ms TCX */
				/* decode gain offset: -10.5, -7.5, -4.5, -1.5, +1.5, +4.5, +7.5, +10.5 dB */
				for( i = 0; i < 16; i++ ) {
					index = *prm_other++;
					tmp = ( 3.0f * ( (float)index ) ) - 10.5f;
					if( bfi ) {
						tmp = 0.0;
					}
					if( bad_frame[k + 1] && ( i < 8 ) ) {
						tmp = 0.0;
					}
					if( bad_frame[k + 2] && ( i >= 8 ) ) {
						tmp = 0.0;
					}
					gain_hf_other[i] = tmp + gain4[( i >> 2 )];
					gain_sum_other += gain_hf_other[i];
				}
			}
			else if( mode == 2 ) { /* 40+5-ms TCX */
				/* decode gain offset: -4.5, -1.5, +1.5, +4.5 dB */
				for( i = 0; i < 8; i++ ) {
					index = *prm_other++;
					tmp = ( 3.0f * ( (float)index ) ) - 4.5f;
					if( bfi || bad_frame[k + 1] ) {
						tmp = 0.0;
					}
					gain_hf_other[i] = tmp + gain4[( i >> 1 )];
					gain_sum_other += gain_hf_other[i];
				}
			}
			else { /* ACELP or 20+2.5-ms TCX */
				for( i = 0; i < 4; i++ ) {
					gain_hf_other[i] = gain4[i];
					gain_sum_other += gain_hf_other[i];
				}
			}
			thr_lo = -20.0f * (float)nsf;
			thr_hi = 20.0f * (float)nsf;
			gdiff = gain_sum - gain_sum_other;
			/* decode ISFs and convert ISFs to cosine domain */
			d_isf_hf( prm_in_other, isfnew_other, st->past_q_isf_hf_other, bfi, st->mean_isf_hf, st->dico1_isf_hf );
			if( gdiff < thr_lo ) {
				for( i = 0; i < MHF; i++ ) {
					isfnew[i] = isfnew_other[i];
				}
			}
			else if( gdiff > thr_hi ) {
			}
			else {
				for( i = 0; i < MHF; i++ ) {
					isfnew[i] = 0.5f * ( isfnew[i] + isfnew_other[i] );
				}
			}
			isf2isp( isfnew, ispnew, MHF );
			int_lpc_np1( st->ispold_hf, ispnew, Aq, nsf, MHF );
			mvr2r( ispnew, st->ispold_hf, MHF );
			/* compute gain correction factor to match amplitude at 6.4kHz */
			gain = match_gain_6k4( &AqLF[( ( k * 4 ) + nsf ) * ( M + 1 )], &Aq[nsf * ( MHF + 1 )] );
			/* interpolate per subframe */
			int_gain( st->old_gain, gain, gain_tmp, nsf );
			for( i = 0; i < nsf; i++ ) {
				gain_hf[i] += gain_tmp[i];
				gain_hf_other[i] += gain_tmp[i];
			}
			st->old_gain = gain;
		}
		for( i = 0; i < MHF + 1; i++ ) {
			mem_lpc_hf[i] = Aq[nsf * ( MHF + 1 ) + i];
		}
		/* HF synthesis and smooting */
		p_Aq = Aq;
		for( sf = 0; sf < nsf; sf++ ) {
			i_subfr = ( ( k * 4 ) + sf ) * L_SUBFR;
			gain = (float)pow( 10.0f, gain_hf[sf] / 20.0f );
			if( mono_dec_stereo == 1 ) {
				gain += (float)pow( 10.0f, gain_hf_other[sf] / 20.0f );
				gain *= 0.5;
			}
			if( *ramp_state < 64 ) {
				gain *= gain_hf_ramp[*ramp_state];
				( *ramp_state )++;
			}
			*mem_gain_hf = gain;
			for( i = 0; i < L_SUBFR; i++ ) {
				HF[i] = gain * exc[i + i_subfr];
			}
			soft_exc_hf( HF, &( st->lp_amp ) );
			/* filtering through HF LP filter */
			E_UTIL_synthesisPlus( p_Aq, MHF, HF, HF, L_SUBFR, st->mem_syn_hf, 1 );
			/* smooth energy of HF signal on subframe (HF[]) */
			smooth_ener_hf( HF, &( st->threshold ) );
			for( i = 0; i < L_SUBFR; i++ ) {
				synth_hf[i + i_subfr] = HF[i];
			}
			p_Aq += ( MHF + 1 );
		}
		if( mode == 2 ) {
			k++;
		}
		else if( mode == 3 ) {
			k += 3;
		}
	}
	return;
}
