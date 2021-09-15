#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
static void try_tcx(
    int    k,
    int    mode, /* 1=TCX20, 2=TCX40 3=TCX80 */
    float *snr,
    float  A[],
    float  wsp[],
    int    mod[],
    short  coding_mod[],
    float  isp[],
    float  isp_q[],
    float  AqLF[],
    float  speech[],
    float  mem_w0[],
    float  mem_xnq[],
    float  mem_wsyn[],
    float  old_exc[],
    float  mem_syn[],
    float  wovlp[],
    int    ovlp_size[],
    float  past_isfq[],
    int    pit_adj,
    int    nbits,
    int    nprm_tcx,
    int    prm[] );
void coder_lf_b( int  codec_mode,   /* (i) : AMR-WB+ mode (see cnst.h)             */
    float             speech[],     /* (i) : speech vector [-M..L_FRAME_PLUS+L_NEXT]    */
    float             synth[],      /* (o) : synthesis vector [-M..L_FRAME_PLUS]        */
    int               mod[],        /* (o) : mode for each 20ms frame (mode[4]     */
    float             AqLF[],       /* (o) : quantized coefficients (AdLF[16])     */
    float             window[],     /* (i) : window for LPC analysis               */
    int               param[],      /* (o) : parameters (NB_DIV*NPRM_DIV)          */
    float             ol_gain[],    /* (o) : open-loop LTP gain                    */
    int               ave_T_out[],  /* (o) : average LTP lag                       */
    float             ave_p_out[],  /* (o) : average LTP gain                      */
    short             coding_mod[], /* (i) : selected mode for each 20ms           */
    int               pit_adj,      /* (i) : indicate pitch adjustment             */
    Coder_State_Plus *st            /* i/o : coder memory state                  */
)
{
	/* LPC coefficients */
	float  r[M + 1]; /* Autocorrelations of windowed speech  */
	float  A[( NB_SUBFR + 1 ) * ( M + 1 )];
	float  Ap[M + 1];
	float  ispnew[M]; /* LSPs at 4nd subframe                 */
	float  isfnew[M];
	float  isp[( NB_DIV + 1 ) * M];
	float  isp_q[( NB_DIV + 1 ) * M];
	float  past_isfq[3 * M]; /* past isf quantizer */
	float  mem_w0[NB_DIV + 1], mem_wsyn[NB_DIV + 1];
	float  mem_xnq[NB_DIV + 1];
	int    ovlp_size[NB_DIV + 1];
	float *wsp;
	float  mem_syn[( NB_DIV + 1 ) * M];
	float  wovlp[( NB_DIV + 1 ) * 128];
	/* Scalars */
	int    i, j, k, i2, i1, nbits, *prm;
	float  snr, snr1, snr2;
	float  ener, cor_max, t0;
	float *p, *p1;
	float  norm_corr[4], norm_corr2[4];
	int    T_op[NB_DIV], T_op2[NB_DIV];
	int    T_out[4]; /* integer pitch-lag */
	float  p_out[4];
	float  LTPGain[2];
	int    PIT_MIN; /* Minimum pitch lag with resolution 1/4      */
	int    PIT_MAX; /* Maximum pitch lag                          */
	if( pit_adj == 0 ) {
		PIT_MIN = PIT_MIN_12k8;
		PIT_MAX = PIT_MAX_12k8;
	}
	else {
		i = ( ( ( pit_adj * PIT_MIN_12k8 ) + ( FSCALE_DENOM / 2 ) ) / FSCALE_DENOM ) - PIT_MIN_12k8;
		PIT_MIN = PIT_MIN_12k8 + i;
		PIT_MAX = PIT_MAX_12k8 + ( 6 * i );
	}
	/* number of bits per frame (80 ms) */
	nbits = NBITS_CORE[codec_mode];
	/* remove bits for mode */
	nbits -= NBITS_MODE;
	/*---------------------------------------------------------------*
	*  Perform LP analysis four times (every 20 ms)                 *
	*  - autocorrelation + lag windowing                            *
	*  - Levinson-Durbin algorithm to find a[]                      *
	*  - convert a[] to isp[]                                       *
	*  - interpol isp[]                                             *
	*---------------------------------------------------------------*/
	/* read old isp for LPC interpolation */
	mvr2r( st->ispold, isp, M );
	for( i = 0; i < NB_DIV; i++ ) {
		if( pit_adj <= FSCALE_DENOM ) {
			/* Autocorrelations of signal at 12.8kHz */
			E_UTIL_autocorrPlus( &speech[( i * L_DIV ) + L_SUBFR], r, M, L_WINDOW, window );
		}
		else {
			/* Autocorrelations of signal at 12.8kHz */
			E_UTIL_autocorrPlus( &speech[( i * L_DIV ) + ( L_SUBFR / 2 )], r, M, L_WINDOW_HIGH_RATE, window );
		}
		/* Lag windowing    */
		lag_wind( r, M );
		/* Levinson Durbin  */
		E_LPC_lev_dur( Ap, r, M );
		/* From A(z) to ISP */
		E_LPC_a_isp_conversion( Ap, ispnew, st->ispold, M );
		for( j = 0; j < M; j++ ) {
			st->stClass->ApBuf[i * M + j] = Ap[j];
		}
		mvr2r( ispnew, &isp[( i + 1 ) * M], M );
		/* A(z) interpolation every 20 ms (used for weighted speech) */
		int_lpc_np1( st->ispold, ispnew, &A[i * 4 * ( M + 1 )], 4, M );
		/* update ispold[] for the next LPC analysis */
		mvr2r( ispnew, st->ispold, M );
	}
	mvr2r( synth - M, mem_syn, M );
	wsp = synth; /* wsp[] use synth[] as tmp buffer */
	{
		float  old_d_wsp[( PIT_MAX_MAX / OPL_DECIM ) + L_DIV]; /* Weighted speech vector */
		float *d_wsp;
		d_wsp = old_d_wsp + PIT_MAX_MAX / OPL_DECIM;
		mvr2r( st->old_d_wsp, old_d_wsp, PIT_MAX_MAX / OPL_DECIM );
		/* Calculate open-loop LTP parameters */
		for( i = 0; i < NB_DIV; i++ ) {
			/* weighted speech for SNR */
			find_wsp( &A[i * ( NB_SUBFR / 4 ) * ( M + 1 )], &speech[i * L_DIV], &wsp[i * L_DIV], &( st->mem_wsp ), L_DIV );
			mvr2r( &wsp[i * L_DIV], d_wsp, L_DIV );
			E_GAIN_lp_decim2( d_wsp, L_DIV, st->mem_lp_decim2 );
			/* Find open loop pitch lag for first 1/2 frame */
			T_op[i] = E_GAIN_open_loop_search( d_wsp, ( PIT_MIN / OPL_DECIM ) + 1, PIT_MAX / OPL_DECIM, ( 2 * L_SUBFR ) / OPL_DECIM, st->old_T0_med, &( st->ol_gain ), st->hp_ol_ltp_mem, st->hp_old_wsp, (unsigned char)st->ol_wght_flg );
			LTPGain[1] = st->ol_gain;
			if( st->ol_gain > 0.6 ) {
				st->old_T0_med = E_GAIN_olag_median( T_op[i], st->old_ol_lag );
				st->ada_w = 1.0;
			}
			else {
				st->ada_w = st->ada_w * 0.9f;
			}
			if( st->ada_w < 0.8 ) {
				st->ol_wght_flg = 0;
			}
			else {
				st->ol_wght_flg = 1;
			}
			/* compute max */
			cor_max = 0.0f;
			p = &d_wsp[0];
			p1 = d_wsp - T_op[i];
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++ ) {
				cor_max += *p++ * *p1++;
			}
			/* compute energy */
			t0 = 0.01f;
			p = d_wsp - T_op[i];
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++, p++ ) {
				t0 += *p * *p;
			}
			t0 = (float)( 1.0 / sqrt( t0 ) ); /* 1/sqrt(energy)    */
			norm_corr[i] = cor_max * t0;      /* max/sqrt(energy)  */
			/* normalized corr (0..1) */
			ener = 0.01f;
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++ ) {
				ener += d_wsp[j] * d_wsp[j];
			}
			ener = (float)( 1.0 / sqrt( ener ) ); /* 1/sqrt(energy)    */
			norm_corr[i] *= ener;
			/* Find open loop pitch lag for second 1/2 frame */
			T_op2[i] = E_GAIN_open_loop_search( d_wsp + ( ( 2 * L_SUBFR ) / OPL_DECIM ), ( PIT_MIN / OPL_DECIM ) + 1, PIT_MAX / OPL_DECIM, ( 2 * L_SUBFR ) / OPL_DECIM, st->old_T0_med, &st->ol_gain, st->hp_ol_ltp_mem, st->hp_old_wsp, (unsigned char)st->ol_wght_flg );
			LTPGain[0] = st->ol_gain;
			if( st->ol_gain > 0.6 ) {
				st->old_T0_med = E_GAIN_olag_median( T_op2[i], st->old_ol_lag );
				st->ada_w = 1.0;
			}
			else {
				st->ada_w = st->ada_w * 0.9f;
			}
			if( st->ada_w < 0.8 ) {
				st->ol_wght_flg = 0;
			}
			else {
				st->ol_wght_flg = 1;
			}
			/* compute max */
			cor_max = 0.0f;
			p = d_wsp + ( 2 * L_SUBFR ) / OPL_DECIM;
			p1 = d_wsp + ( ( 2 * L_SUBFR ) / OPL_DECIM ) - T_op2[i];
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++ ) {
				cor_max += *p++ * *p1++;
			}
			/* compute energy */
			t0 = 0.01f;
			p = d_wsp + ( ( 2 * L_SUBFR ) / OPL_DECIM ) - T_op2[i];
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++, p++ ) {
				t0 += *p * *p;
			}
			t0 = (float)( 1.0 / sqrt( t0 ) ); /* 1/sqrt(energy)    */
			norm_corr2[i] = cor_max * t0;     /* max/sqrt(energy)  */
			/* normalized corr (0..1) */
			ener = 0.01f;
			for( j = 0; j < ( 2 * L_SUBFR ) / OPL_DECIM; j++ ) {
				ener += d_wsp[( ( 2 * L_SUBFR ) / OPL_DECIM ) + j] * d_wsp[( ( 2 * L_SUBFR ) / OPL_DECIM ) + j];
			}
			ener = (float)( 1.0 / sqrt( ener ) ); /* 1/sqrt(energy)    */
			norm_corr2[i] *= ener;
			ol_gain[i] = st->ol_gain;
			mvr2r( &old_d_wsp[L_DIV / OPL_DECIM], old_d_wsp, PIT_MAX_MAX / OPL_DECIM );
			st->stClass->LTPGain[( 2 * i + 2 )] = LTPGain[1];
			st->stClass->LTPGain[( 2 * i + 2 ) + 1] = LTPGain[0];
			st->stClass->LTPLag[( 2 * i ) + 2] = T_op[i];
			st->stClass->LTPLag[( 2 * i + 2 ) + 1] = T_op2[i];
			st->stClass->NormCorr[( 2 * i + 2 )] = norm_corr[i];
			st->stClass->NormCorr[( 2 * i + 2 ) + 1] = norm_corr2[i];
		}
		mvr2r( old_d_wsp, st->old_d_wsp, PIT_MAX_MAX / OPL_DECIM ); /* d_wsp already shifted */
	}
	/*---------------------------------------------------------------*
  *  Call ACELP and TCX codec                                     *
  *---------------------------------------------------------------*/
	ovlp_size[0] = st->old_ovlp_size;
	/*Classification refinement */
	classifyExcitationRef( st->stClass, isp, coding_mod );
	if( st->SwitchFlagPlusToWB && coding_mod[0] != 0 ) {
		coding_mod[0] = 0;
		st->stClass->NbOfAcelps++;
		st->SwitchFlagPlusToWB = 0;
	}
	ovlp_size[1] = st->old_ovlp_size;
	ovlp_size[2] = st->old_ovlp_size;
	ovlp_size[3] = st->old_ovlp_size;
	ovlp_size[4] = st->old_ovlp_size;
	mem_w0[0] = st->old_mem_w0;
	mem_xnq[0] = st->old_mem_xnq;
	mem_wsyn[0] = st->old_mem_wsyn;
	mvr2r( st->old_wovlp, wovlp, 128 );
	mvr2r( st->past_isfq, &past_isfq[0], M );
	mvr2r( st->ispold_q, isp_q, M );
	snr2 = 0.0;
	for( i1 = 0; i1 < 2; i1++ ) {
		mvr2r( &past_isfq[i1 * M], &past_isfq[( i1 + 1 ) * M], M );
		snr1 = 0.0;
		for( i2 = 0; i2 < 2; i2++ ) {
			k = ( i1 * 2 ) + i2;
			if( coding_mod[k] == 0 || ( coding_mod[k] == 1 && st->stClass->NbOfAcelps != 0 )
			    || ( st->stClass->NoMtcx[i1] != 0 && st->stClass->NbOfAcelps == 0 ) ) {
				/* set pointer to parameters */
				prm = param + ( k * NPRM_DIV );
				/*---------------------------------------------------------------*
			  *  Quantize ISF parameters (46 bits) every 20 ms                *
			  *---------------------------------------------------------------*/
				/* Convert isps to frequency domain 0..6400 */
				E_LPC_isp_isf_conversion( &isp[( k + 1 ) * M], isfnew, M );
				/* quantize 1st and 2nd LPC with 46 bits */
				qpisf_2s_46b( isfnew, isfnew, &past_isfq[( i1 + 1 ) * M], prm, 4 );
				prm += NPRM_LPC;
				/* Convert isfs to the cosine domain */
				isf2isp( isfnew, &isp_q[( k + 1 ) * M], M );
				/* interpol quantized lpc */
				int_lpc_np1( &isp_q[k * M], &isp_q[( k + 1 ) * M], &AqLF[k * 4 * ( M + 1 )], 4, M );
				/*---------------------------------------------------------------*
			  *  Call ACELP 4 subfr x 5ms = 20 ms frame                       *
			  *---------------------------------------------------------------*/
				if( coding_mod[k] == 0 ) {
					mem_w0[k + 1] = mem_w0[k];
					mem_xnq[k + 1] = mem_xnq[k];
					ovlp_size[k + 1] = 0;
					{
						float old_exc[PIT_MAX_MAX + L_INTERPOL + L_DIV + 1]; /* Excitation vector */
						float old_syn[M + L_DIV];
						float buf[L_DIV];
						mvr2r( st->old_exc, old_exc, PIT_MAX_MAX + L_INTERPOL );
						mvr2r( &mem_syn[k * M], old_syn, M );
						coder_acelp( &A[k * ( NB_SUBFR / 4 ) * ( M + 1 )],
						    &AqLF[k * ( NB_SUBFR / 4 ) * ( M + 1 )],
						    &speech[k * L_DIV],
						    &mem_w0[k + 1],
						    &mem_xnq[k + 1],
						    old_syn + M,
						    old_exc + PIT_MAX_MAX + L_INTERPOL,
						    &wovlp[( k + 1 ) * 128],
						    L_DIV,
						    codec_mode,
						    norm_corr[k],
						    norm_corr2[k],
						    T_op[k],
						    T_op2[k],
						    T_out,
						    p_out,
						    st->mem_gain_code,
						    pit_adj,
						    prm );
						mvr2r( &old_exc[L_DIV], st->old_exc, PIT_MAX_MAX + L_INTERPOL );
						mvr2r( old_syn + L_DIV, &mem_syn[( k + 1 ) * M], M );
						/* average integer pitch-lag for high band coding */
						ave_T_out[k] = (int)( ( T_op[k] + T_op2[k] ) );
						ave_p_out[k] = ( p_out[0] + p_out[1] + p_out[2] + p_out[3] ) / 4.0f;
						mem_wsyn[k + 1] = mem_wsyn[k];
						find_wsp( &A[k * ( NB_SUBFR / 4 ) * ( M + 1 )], old_syn + M, buf, &mem_wsyn[k + 1], L_DIV );
					}
					mod[k] = 0;
					coding_mod[k] = 0;
				} /* end ACELP */
				  /*--------------------------------------------------*
				* Call 20MS TCX coder and find segmental SNR       *
			  *--------------------------------------------------*/
				if( coding_mod[k] != 0 ) {
					snr = -10000.0f;
					try_tcx( k, 1, &snr, A, wsp, mod, coding_mod, isp, isp_q, AqLF, speech, mem_w0, mem_xnq, mem_wsyn, st->old_exc, mem_syn, wovlp, ovlp_size, past_isfq, pit_adj, ( ( nbits / 4 ) - NBITS_LPC ), NPRM_TCX20, prm );
					snr1 += 0.5f * snr;
				} /* end of TCX_20MS */
			}     /* end of i2 */
		}         /* end of coding mode ACELP or TCX_20MS */
		if( coding_mod[i1 * 2] != 0 && coding_mod[i1 * 2 + 1] != 0 ) {
			if( st->stClass->NbOfAcelps == 0 ) {
				snr1 = -10000.0; /* TCX20 off*/
			}
			k = ( i1 * 2 );
			/* set pointer to parameters */
			prm = param + ( k * NPRM_DIV );
			try_tcx( k, 2, &snr1, A, wsp, mod, coding_mod, isp, isp_q, AqLF, speech, mem_w0, mem_xnq, mem_wsyn, st->old_exc, mem_syn, wovlp, ovlp_size, past_isfq, pit_adj, ( ( nbits / 2 ) - NBITS_LPC ), NPRM_LPC + NPRM_TCX40, prm );
			snr2 += 0.5f * snr1;
		} /* end of coding mode TCX_40MS */
	}     /* end of i1 */
	if( coding_mod[0] != 0 && coding_mod[1] != 0 && coding_mod[2] != 0 && coding_mod[3] != 0 && st->stClass->NoMtcx[0] == 0 && st->stClass->NoMtcx[1] == 0 ) {
		k = 0;
		/* set pointer to parameters */
		prm = param + ( k * NPRM_DIV );
		try_tcx( k, 3, &snr2, A, wsp, mod, coding_mod, isp, isp_q, AqLF, speech, mem_w0, mem_xnq, mem_wsyn, st->old_exc, mem_syn, wovlp, ovlp_size, past_isfq, pit_adj, ( nbits - NBITS_LPC ), NPRM_LPC + NPRM_TCX80, prm );
	} /* end of coding mode 3 */
	/*Mode buffer update & reset MTcx counter*/
	for( i = 0; i < 4; i++ ) {
		st->stClass->prevModes[i] = coding_mod[i];
	}
	if( st->stClass->NoMtcx[0] != 0 ) {
		st->stClass->NoMtcx[0] = 0;
	}
	if( st->stClass->NoMtcx[1] != 0 ) {
		st->stClass->NoMtcx[1] = 0;
	}
	/*--------------------------------------------------*
	 * Update filter memory.                            *
	 *--------------------------------------------------*/
	st->old_ovlp_size = ovlp_size[NB_DIV];
	st->old_mem_w0 = mem_w0[NB_DIV];
	st->old_mem_xnq = mem_xnq[NB_DIV];
	st->old_mem_wsyn = mem_wsyn[NB_DIV];
	mvr2r( &wovlp[NB_DIV * 128], st->old_wovlp, 128 );
	mvr2r( &past_isfq[2 * M], st->past_isfq, M );
	mvr2r( &isp_q[NB_DIV * M], st->ispold_q, M );
	mvr2r( &mem_syn[NB_DIV * M], synth + L_FRAME_PLUS - M, M );
	/*--------------------------------------------------*
	 * Update exc for next frame.                       *
	 *--------------------------------------------------*/
	return;
}
static void try_tcx(
    int    k,
    int    mode, /* 1=TCX20, 2=TCX40 3=TCX80 */
    float *snr,
    float  A[],
    float  wsp[],
    int    mod[],
    short  coding_mod[],
    float  isp[],
    float  isp_q[],
    float  AqLF[],
    float  speech[],
    float  mem_w0[],
    float  mem_xnq[],
    float  mem_wsyn[],
    float  old_exc[],
    float  mem_syn[],
    float  wovlp[],
    int    ovlp_size[],
    float  past_isfq[],
    int    pit_adj,
    int    nbits,
    int    nprm_tcx,
    int    prm[] )
{
	float Aq[( NB_SUBFR + 1 ) * ( M + 1 )], *p_Aq;
	int   prm_tcx_buf[NPRM_LPC + NPRM_TCX80], *prm_tcx;
	float synth_tcx[M + L_TCX];
	float exc_tcx[L_TCX];
	float mem_w0_tcx, mem_xnq_tcx, mem_wsyn_tcx;
	float wovlp_tcx[128];
	float past_isfq_tcx[M];
	float ispnew_q[M]; /* LSPs at 4nd subframe                 */
	float isfnew[M];
	/* Scalars */
	int   i, ndiv;
	float tmp;
	int   PIT_MAX;
	if( pit_adj == 0 ) {
		PIT_MAX = PIT_MAX_12k8;
	}
	else {
		i = ( ( ( pit_adj * PIT_MIN_12k8 ) + ( FSCALE_DENOM / 2 ) ) / FSCALE_DENOM ) - PIT_MIN_12k8;
		PIT_MAX = PIT_MAX_12k8 + ( 6 * i );
	}
	prm_tcx = prm_tcx_buf;
	ndiv = mode;
	if( ndiv == 3 ) {
		ndiv++; /* 4 divisions in mode 3 (TCX80) */
	}
	if( mode > 1 ) {
		/*---------------------------------------------------------------*
    *  Quantize ISF parameters (46 bits) every 40/80 ms             *
    *---------------------------------------------------------------*/
		/* Convert isps to frequency domain 0..6400 */
		E_LPC_isp_isf_conversion( &isp[( k + ndiv ) * M], isfnew, M );
		mvr2r( &past_isfq[( k / 2 ) * M], past_isfq_tcx, M );
		/* quantize 1st and 2nd LPC with 46 bits */
		qpisf_2s_46b( isfnew, isfnew, past_isfq_tcx, prm_tcx, 4 );
		prm_tcx += NPRM_LPC;
		/* Convert isfs to the cosine domain */
		isf2isp( isfnew, ispnew_q, M );
		/* interpol quantized lpc */
		int_lpc_np1( &isp_q[k * M], ispnew_q, Aq, ndiv * 4, M );
		p_Aq = Aq;
	}
	else {
		p_Aq = &AqLF[k * 4 * ( M + 1 )];
	}
	/*--------------------------------------------------------*
  * Call 20/40/80MS TCX coder and find segmental SNR       *
  *--------------------------------------------------------*/
	mvr2r( &mem_syn[k * M], synth_tcx, M );
	mem_w0_tcx = mem_w0[k];
	mem_xnq_tcx = mem_xnq[k];
	mvr2r( &wovlp[k * 128], wovlp_tcx, 128 );
	coder_tcx( p_Aq,
	    &speech[k * L_DIV],
	    &mem_w0_tcx,
	    &mem_xnq_tcx,
	    &synth_tcx[M],
	    exc_tcx,
	    wovlp_tcx,
	    ovlp_size[k],
	    ndiv * L_DIV,
	    nbits,
	    prm_tcx );
	mem_wsyn_tcx = mem_wsyn[k];
	{
		float buf[L_FRAME_PLUS];
		find_wsp( &A[k * 4 * ( M + 1 )], &synth_tcx[M], buf, &mem_wsyn_tcx, ndiv * L_DIV );
		tmp = segsnr( &wsp[k * L_DIV], buf, (short)( ndiv * L_DIV ), (short)L_SUBFR );
	}
	/*--------------------------------------------------------*
  * Save tcx parameters if tcx segmental SNR is better     *
  *--------------------------------------------------------*/
	if( tmp > *snr ) {
		*snr = tmp;
		for( i = 0; i < ndiv; i++ ) {
			mod[k + i] = mode;
			coding_mod[k + i] = mode;
		}
		mem_w0[k + ndiv] = mem_w0_tcx;
		mem_xnq[k + ndiv] = mem_xnq_tcx;
		mem_wsyn[k + ndiv] = mem_wsyn_tcx;
		ovlp_size[k + ndiv] = ( 32 * ndiv );
		if( mode > 1 ) {
			mvr2r( ispnew_q, &isp_q[( k + ndiv ) * M], M );
			mvr2r( past_isfq_tcx, &past_isfq[( ( k + ndiv ) / 2 ) * M], M );
			/* lpc coefficient needed for HF extension */
			mvr2r( Aq, &AqLF[k * 4 * ( M + 1 )], ( ( ndiv * 4 ) + 1 ) * ( M + 1 ) );
		}
		if( ndiv == 1 ) {
			mvr2r( &old_exc[L_DIV], old_exc, ( PIT_MAX_MAX + L_INTERPOL ) - L_DIV );
			mvr2r( exc_tcx, old_exc + PIT_MAX_MAX + L_INTERPOL - L_DIV, L_DIV );
		}
		else {
			mvr2r( exc_tcx + ( ndiv * L_DIV ) - ( PIT_MAX_MAX + L_INTERPOL ), old_exc, PIT_MAX_MAX + L_INTERPOL );
		}

		mvr2r( wovlp_tcx, &wovlp[( k + ndiv ) * 128], 128 );
		mvr2r( &synth_tcx[ndiv * L_DIV], &mem_syn[( k + ndiv ) * M], M );
		mvi2i( prm_tcx_buf, prm, nprm_tcx );
	}
	return;
}
