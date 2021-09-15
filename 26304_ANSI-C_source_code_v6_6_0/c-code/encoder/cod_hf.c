#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*----------- CODE AND SYNTHESIZE UPPER-BAND SIGNAL (synth_hf[])  ----------
    bitstream pointer = ptr
  coding modes (coder type=ACELP,TCX and frame length=20,40,80 ms) for upper-band are copied from lower-band
  HF coding mode depends on bit allocation to HFs (nbits_hf) :
  nbits_hf = 4*16=64
       the lower-band excitation is randomized, folded, and shaped first in time-domain
       with subframe gains then in frequency domain with a LP filter (plus energy smoothing)
  nbits_hf = 4*16+640=704
       the lower-band is quantized by ACELP/TCX model, using the same mode as lower band.
   */
/*------------------------------------------------------*
  * HF analysis, quantization & synthesis                *
  *------------------------------------------------------*/
void init_coder_hf( Coder_StState *st )
{
	int i;
	/* Static vectors to zero */
	set_zero( st->mem_hf1, MHF );
	set_zero( st->mem_hf2, MHF );
	set_zero( st->mem_hf3, MHF );
	set_zero( st->past_q_isf_hf, MHF );
	st->old_gain = 1.0;
	/* isp initialization */
	for( i = 0; i < MHF - 1; i++ )
		st->ispold_hf[i] = (float)cos( 3.141592654 * (float)( i + 1 ) / (float)MHF );
	st->ispold_hf[MHF - 1] = 0.045f;
	mvr2r( st->ispold_hf, st->ispold_q_hf, MHF );
	return;
}
void coder_hf(
    int            mod[],       /* (i) : mode for each 20ms frame (mode[4]     */
    float          AqLF[],      /* (i) : Q coeff in lower band (AdLF[16])      */
    float          speech[],    /* (i) : speech vector [-M..L_FRAME_PLUS]           */
    float          speech_hf[], /* (i) : HF speech vec [-MHF..L_FRAME_PLUS+L_NEXT]  */
    float          synth_hf[],  /* (o) : HF synthesis  [0..L_FRAME_PLUS]            */
    float          window[],    /* (i) : window for LPC analysis               */
    int            param[],     /* (o) : parameters (NB_DIV*NPRM_HF_DIV)       */
    int            fscale,
    Coder_StState *st ) /* i/o : coder memory state                    */
{
	/* LPC coefficients */
	float  r[MHF + 1]; /* Autocorrelations of windowed speech  */
	float  A[( NB_SUBFR + 1 ) * ( MHF + 1 )];
	float  Aq[( NB_SUBFR + 1 ) * ( MHF + 1 )];
	float *p_A, *p_Aq, Ap[MHF + 1];
	float  ispnew[NB_DIV * MHF]; /* LSPs at 4nd subframe                 */
	float  ispnew_q[MHF];        /* LSPs at 4nd subframe                 */
	float  isfnew[MHF];
	float  gain_hf[NB_SUBFR], gain4[4];
	float  HF[L_SUBFR];
	float  buf[MHF + L_SUBFR];
	/* Scalars */
	int   i, j, k, i_subfr, sf, nsf, index, mode, *prm;
	float tmp, gain, ener;
	/*---------------------------------------------------------------*
  *  Perform HF LP analysis four times (every 20 ms)              *
  *  - autocorrelation + lag windowing                            *
  *  - Levinson-Durbin algorithm to find a[]                      *
  *  - convert a[] to isp[]                                       *
  *  - interpol isp[]                                             *
  *---------------------------------------------------------------*/
	for( k = 0; k < NB_DIV; k++ ) {
		if( fscale <= FSCALE_DENOM ) {
			/* Autocorrelations of HF signal folded into 12.8kHz */
			E_UTIL_autocorrPlus( &speech_hf[( k * L_DIV ) + L_SUBFR], r, MHF, L_WINDOW, window );
		}
		else {
			/* Autocorrelations of HF signal folded into 12.8kHz */
			E_UTIL_autocorrPlus( &speech_hf[( k * L_DIV ) + ( L_SUBFR / 2 )], r, MHF, L_WINDOW_HIGH_RATE, window );
		}
		lag_wind( r, MHF );
		E_LPC_lev_dur( Ap, r, MHF );
		E_LPC_a_isp_conversion( Ap, &ispnew[k * MHF], st->ispold_hf, MHF );
		int_lpc_np1( st->ispold_hf, &ispnew[k * MHF], &A[k * 4 * ( MHF + 1 )], 4, MHF );
		mvr2r( &ispnew[k * MHF], st->ispold_hf, MHF );
	}
	p_A = A;
	for( k = 0; k < NB_DIV; k++ ) {
		mode = mod[k];
		if( mode == 0 ) {
			nsf = 4;
		}
		else {
			nsf = 2 << mode;
		}
		/* set pointer to parameters */
		prm = param + ( k * NPRM_BWE_DIV );
		/* Convert isps to frequency domain 0..6400 */
		E_LPC_isp_isf_conversion( &ispnew[( k + ( nsf / 4 ) - 1 ) * MHF], isfnew, MHF );
		/* quantize isf */
		q_isf_hf( isfnew, isfnew, st->past_q_isf_hf, prm, st->mean_isf_hf, st->dico1_isf_hf );
		prm += 2;
		/* Convert isfs to the cosine domain */
		isf2isp( isfnew, ispnew_q, MHF );
		/* interpol quantized isp for local synthesis */
		int_lpc_np1( st->ispold_q_hf, ispnew_q, Aq, nsf, MHF );
		mvr2r( ispnew_q, st->ispold_q_hf, MHF );
		/* gain factor to match amplitude at 6.4kHz */
		gain = match_gain_6k4( &AqLF[( ( k * 4 ) + nsf ) * ( M + 1 )], &Aq[nsf * ( MHF + 1 )] );
		int_gain( st->old_gain, gain, gain_hf, nsf );
		st->old_gain = gain;
		/*------------------------------------------------------*
    * find gain in weighted domain                         *
    *------------------------------------------------------*/
		p_Aq = Aq;
		for( sf = 0; sf < nsf; sf++ ) {
			i_subfr = ( ( k * 4 ) + sf ) * L_SUBFR;
			E_UTIL_residu( &AqLF[( ( k * 4 ) + sf ) * ( M + 1 )], &speech[i_subfr], buf + MHF, L_SUBFR );
			mvr2r( st->mem_hf1, buf, MHF );
			E_UTIL_synthesisPlus( p_Aq, MHF, buf + MHF, buf + MHF, L_SUBFR, st->mem_hf1, 1 );
			E_LPC_a_weight( p_A, Ap, GAMMA_HF, MHF );
			E_UTIL_residuPlus( p_A, MHF, buf + MHF, HF, L_SUBFR );
			E_UTIL_synthesisPlus( Ap, MHF, HF, HF, L_SUBFR, st->mem_hf2, 1 );
			E_UTIL_residuPlus( p_A, MHF, &speech_hf[i_subfr], buf, L_SUBFR );
			E_UTIL_synthesisPlus( Ap, MHF, buf, buf, L_SUBFR, st->mem_hf3, 1 );
			ener = 0.001f;
			for( i = 0; i < L_SUBFR; i++ ) {
				ener += buf[i] * buf[i];
			}
			tmp = 0.001f;
			for( i = 0; i < L_SUBFR; i++ ) {
				tmp += HF[i] * HF[i];
			}
			gain = 10.0f * (float)log10( ener / tmp );
			gain_hf[sf] = gain - gain_hf[sf];
			p_A += ( MHF + 1 );
			p_Aq += ( MHF + 1 );
		}
		/*------------------------------------------------------*
    * HF gains quantization                                *
    *------------------------------------------------------*/
		if( mode == 3 ) /* TCX 80 */
		{
			for( i = 0; i < 4; i++ ) {
				tmp = 0.0;
				for( j = 0; j < 4; j++ ) {
					tmp += gain_hf[( 4 * i ) + j];
				}
				gain4[i] = 0.25f * tmp;
			}
			q_gain_hf( gain4, gain4, prm );
			prm++;
			for( i = 0; i < nsf; i++ ) {
				gain = gain4[i >> 2];
				/* quantize gain at: -10.5, -7.5, -4.5, -1.5, +1.5, +4.5, +7.5, +10.5 dB */
				tmp = gain_hf[i] - gain;
				index = (int)floor( ( ( 10.5 + tmp ) / 3.0 ) + 0.5 );
				if( index > 7 ) {
					index = 7;
				}
				if( index < 0 ) {
					index = 0;
				}
				prm[i] = index;
				gain_hf[i] = gain + ( 3.0f * ( (float)index ) ) - 10.5f;
			}
			prm += nsf;
		}
		else if( mode == 2 ) /* TCX 40 */
		{
			for( i = 0; i < 4; i++ ) {
				tmp = 0.0;
				for( j = 0; j < 2; j++ ) {
					tmp += gain_hf[( 2 * i ) + j];
				}
				gain4[i] = 0.5f * tmp;
			}
			q_gain_hf( gain4, gain4, prm );
			prm++;
			for( i = 0; i < nsf; i++ ) {
				gain = gain4[i >> 1];
				/* quantize gain at: -4.5, -1.5, +1.5, +4.5 dB */
				tmp = gain_hf[i] - gain;
				index = (int)floor( ( ( 4.5 + tmp ) / 3.0 ) + 0.5 );
				if( index > 3 ) {
					index = 3;
				}
				if( index < 0 ) {
					index = 0;
				}
				prm[i] = index;
				gain_hf[i] = gain + ( 3.0f * ( (float)index ) ) - 4.5f;
			}
			prm += nsf;
		}
		else /* ACELP/TCX 20 */
		{
			q_gain_hf( gain_hf, gain_hf, prm );
			prm++;
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
