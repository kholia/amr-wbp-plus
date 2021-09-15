#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*-----------------------------------------------------------------*
 *   Funtion init_coder_wb
 *   ~~~~~~~~~~~~~
 * ->Initialization of variables for the coder section.
 * - initilize pointers to speech buffer
 * - initialize static pointers
 * - set static vectors to zero
 * - compute lag window and LP analysis window
 *
 *-----------------------------------------------------------------*/
void init_coder_lf( Coder_State_Plus *st )
{
	int i;
	/* Static vectors to zero */
	set_zero( st->old_exc, PIT_MAX_MAX + L_INTERPOL );
	set_zero( st->old_d_wsp, PIT_MAX_MAX / OPL_DECIM );
	set_zero( st->mem_lp_decim2, 3 );
	set_zero( st->past_isfq, M );
	set_zero( st->old_wovlp, 128 );
	set_zero( st->hp_old_wsp, L_FRAME_PLUS / OPL_DECIM + ( PIT_MAX_MAX / OPL_DECIM ) );
	set_zero( st->hp_ol_ltp_mem, 3 * 2 + 1 );
	for( i = 0; i < 5; i++ )
		st->old_ol_lag[i] = 40;
	st->old_mem_wsyn = 0.0;
	st->old_mem_w0 = 0.0;
	st->old_mem_xnq = 0.0;
	st->mem_wsp = 0.0;
	st->old_ovlp_size = 0;
	st->old_T0_med = 0;
	st->ol_wght_flg = 0;
	st->ada_w = 0.0;
	/* isf and isp initialization */
	mvr2r( (float *)isf_init, st->isfold, M );
	for( i = 0; i < M - 1; i++ )
		st->ispold[i] = (float)cos( 3.141592654 * (float)( i + 1 ) / (float)M );
	st->ispold[M - 1] = 0.045f;
	mvr2r( st->ispold, st->ispold_q, M );
	return;
}
/*-----------------------------------------------------------------*
 *   Funtion  coder_wb                                             *
 *            ~~~~~~~~                                             *
 *   ->Main coder routine.                                         *
 *                                                                 *
 *-----------------------------------------------------------------*/
void coder_lf(
    int               codec_mode,   /* (i) : AMR-WB+ mode (see cnst.h)             */
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
    Coder_State_Plus *st            /* i/o : coder memory state               */
)
{
	/* LPC coefficients */
	float  r[M + 1]; /* Autocorrelations of windowed speech  */
	float  A[( NB_SUBFR + 1 ) * ( M + 1 )];
	float  Aq[( NB_SUBFR + 1 ) * ( M + 1 )];
	float  ispnew[M];   /* LSPs at 4nd subframe                 */
	float  ispnew_q[M]; /* LSPs at 4nd subframe                 */
	float  isp[( NB_DIV + 1 ) * M];
	float  isp_q[( NB_DIV + 1 ) * M];
	float  isfnew[M];
	float  past_isfq_1[M]; /* past isf quantizer */
	float  past_isfq_2[M]; /* past isf quantizer */
	float  mem_w0[NB_DIV + 1], mem_wsyn[NB_DIV + 1];
	float  mem_xnq[NB_DIV + 1];
	int    ovlp_size[NB_DIV + 1];
	float  Ap[M + 1];
	int    prm_tcx[NPRM_LPC + NPRM_TCX80];
	float  synth_tcx[M + L_TCX];
	float  exc_tcx[L_TCX], mem_w0_tcx, mem_xnq_tcx, mem_wsyn_tcx;
	float  wsp[L_FRAME_PLUS];
	float  wovlp[( NB_DIV + 1 ) * 128];
	float  wovlp_tcx[128];
	float  old_d_wsp[( PIT_MAX_MAX / OPL_DECIM ) + L_DIV];       /* Weighted speech vector */
	float  old_exc[PIT_MAX_MAX + L_INTERPOL + L_FRAME_PLUS + 1]; /* Excitation vector */
	float *d_wsp, *exc;
	/* Scalars */
	int    i, j, k, i2, i1, nbits, *prm;
	float  snr, snr1, snr2;
	float  tmp;
	float  ener, cor_max, t0;
	float *p, *p1;
	float  norm_corr[4], norm_corr2[4];
	int    T_op[NB_DIV], T_op2[NB_DIV];
	int    T_out[4]; /* integer pitch-lag */
	float  p_out[4];
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
	/* Initialize pointers */
	d_wsp = old_d_wsp + PIT_MAX_MAX / OPL_DECIM;
	exc = old_exc + PIT_MAX_MAX + L_INTERPOL;
	/* copy coder memory state into working space (dynamic memory) */
	mvr2r( st->old_d_wsp, old_d_wsp, PIT_MAX_MAX / OPL_DECIM );
	mvr2r( st->old_exc, old_exc, PIT_MAX_MAX + L_INTERPOL );
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
		lag_wind( r, M );                                    /* Lag windowing    */
		E_LPC_lev_dur( Ap, r, M );                           /* Levinson Durbin  */
		E_LPC_a_isp_conversion( Ap, ispnew, st->ispold, M ); /* From A(z) to ISP */
		mvr2r( ispnew, &isp[( i + 1 ) * M], M );
		/* A(z) interpolation every 20 ms (used for weighted speech) */
		int_lpc_np1( st->ispold, ispnew, &A[i * 4 * ( M + 1 )], 4, M );
		/* update ispold[] for the next LPC analysis */
		mvr2r( ispnew, st->ispold, M );
	}
	/*---------------------------------------------------------------*
   * Calculate open-loop LTP parameters                            *
   *---------------------------------------------------------------*/
	for( i = 0; i < NB_DIV; i++ ) {
		/* weighted speech for SNR */
		find_wsp( &A[i * ( NB_SUBFR / 4 ) * ( M + 1 )], &speech[i * L_DIV], &wsp[i * L_DIV], &( st->mem_wsp ), L_DIV );
		mvr2r( &wsp[i * L_DIV], d_wsp, L_DIV );
		E_GAIN_lp_decim2( d_wsp, L_DIV, st->mem_lp_decim2 );
		/* Find open loop pitch lag for first 1/2 frame */
		T_op[i] = E_GAIN_open_loop_search( d_wsp, ( PIT_MIN / OPL_DECIM ) + 1, PIT_MAX / OPL_DECIM, ( 2 * L_SUBFR ) / OPL_DECIM, st->old_T0_med, &( st->ol_gain ), st->hp_ol_ltp_mem, st->hp_old_wsp, (unsigned char)st->ol_wght_flg );
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
	}
	/*---------------------------------------------------------------*
   *  Call ACELP and TCX codec                                     *
   *---------------------------------------------------------------*/
	ovlp_size[0] = st->old_ovlp_size;
	mem_w0[0] = st->old_mem_w0;
	mem_xnq[0] = st->old_mem_xnq;
	mem_wsyn[0] = st->old_mem_wsyn;
	mvr2r( st->old_wovlp, wovlp, 128 );
	mvr2r( st->past_isfq, past_isfq_1, M );
	mvr2r( st->ispold_q, isp_q, M );
	snr2 = 0.0;
	for( i1 = 0; i1 < 2; i1++ ) {
		mvr2r( past_isfq_1, past_isfq_2, M );
		snr1 = 0.0;
		for( i2 = 0; i2 < 2; i2++ ) {
			k = ( i1 * 2 ) + i2;
			/* set pointer to parameters */
			prm = param + ( k * NPRM_DIV );
			/*---------------------------------------------------------------*
       *  Quantize ISF parameters (46 bits) every 20 ms                *
       *---------------------------------------------------------------*/
			/* Convert isps to frequency domain 0..6400 */
			E_LPC_isp_isf_conversion( &isp[( k + 1 ) * M], isfnew, M );
			/* quantize 1st and 2nd LPC with 46 bits */
			qpisf_2s_46b( isfnew, isfnew, past_isfq_1, prm, 4 );
			prm += NPRM_LPC;
			/* Convert isfs to the cosine domain */
			isf2isp( isfnew, &isp_q[( k + 1 ) * M], M );
			/* interpol quantized lpc */
			int_lpc_np1( &isp_q[k * M], &isp_q[( k + 1 ) * M], Aq, 4, M );
			/* lpc coefficient needed for HF extension */
			mvr2r( Aq, &AqLF[k * 4 * ( M + 1 )], 5 * ( M + 1 ) );
			/*---------------------------------------------------------------*
       *  Call ACELP 4 subfr x 5ms = 20 ms frame                       *
       *---------------------------------------------------------------*/
			mem_w0[k + 1] = mem_w0[k];
			mem_xnq[k + 1] = mem_xnq[k];
			ovlp_size[k + 1] = 0;
			coder_acelp( &A[k * ( NB_SUBFR / 4 ) * ( M + 1 )],
			    Aq,
			    &speech[k * L_DIV],
			    &mem_w0[k + 1],
			    &mem_xnq[k + 1],
			    &synth[k * L_DIV],
			    &exc[k * L_DIV],
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
			/* average integer pitch-lag for high band coding */
			ave_T_out[k] = (int)( ( T_op[k] + T_op2[k] ) );
			ave_p_out[k] = ( p_out[0] + p_out[1] + p_out[2] + p_out[3] ) / 4.0f;
			mem_wsyn[k + 1] = mem_wsyn[k];
			{
				float buf[L_FRAME_PLUS];
				find_wsp( &A[k * ( NB_SUBFR / 4 ) * ( M + 1 )], &synth[k * L_DIV], buf, &mem_wsyn[k + 1], L_DIV );
				snr = segsnr( &wsp[k * L_DIV], buf, L_DIV, L_SUBFR );
				if( st->SwitchFlagPlusToWB ) {
					snr = 100000;
					st->SwitchFlagPlusToWB = 0;
				}
			}
			mod[k] = 0;
			coding_mod[k] = 0;
			/*--------------------------------------------------*
       * Call short TCX coder and find segmental SNR       *
       *--------------------------------------------------*/
			mvr2r( &synth[( k * L_DIV ) - M], synth_tcx, M );
			mem_w0_tcx = mem_w0[k];
			mem_xnq_tcx = mem_xnq[k];
			mvr2r( &wovlp[k * 128], wovlp_tcx, 128 );
			coder_tcx( Aq,
			    &speech[k * L_DIV],
			    &mem_w0_tcx,
			    &mem_xnq_tcx,
			    &synth_tcx[M],
			    exc_tcx,
			    wovlp_tcx,
			    ovlp_size[k],
			    L_DIV,
			    ( ( nbits / 4 ) - NBITS_LPC ),
			    prm_tcx );
			mem_wsyn_tcx = mem_wsyn[k];
			{
				float buf[L_FRAME_PLUS];
				find_wsp( &A[k * ( NB_SUBFR / 4 ) * ( M + 1 )], &synth_tcx[M], buf, &mem_wsyn_tcx, L_DIV );
				tmp = segsnr( &wsp[k * L_DIV], buf, L_DIV, L_SUBFR );
			}
			/*--------------------------------------------------------*
       * Save tcx parameters if tcx segmental SNR is better     *
       *--------------------------------------------------------*/
			if( tmp > snr ) {
				snr = tmp;
				mod[k] = 1;
				coding_mod[k] = 1;
				mem_w0[k + 1] = mem_w0_tcx;
				mem_xnq[k + 1] = mem_xnq_tcx;
				mem_wsyn[k + 1] = mem_wsyn_tcx;
				ovlp_size[k + 1] = 32;
				mvr2r( wovlp_tcx, &wovlp[( k + 1 ) * 128], 128 );
				mvr2r( &synth_tcx[M], &synth[k * L_DIV], L_DIV );
				mvr2r( exc_tcx, &exc[k * L_DIV], L_DIV );
				mvi2i( prm_tcx, prm, NPRM_TCX20 );
			}
			snr1 += 0.5f * snr;
		} /* end of coding mode ACELP or TCX_20MS */
		k = ( i1 * 2 );
		/* set pointer to parameters */
		prm = param + ( k * NPRM_DIV );
		/*---------------------------------------------------------------*
     *  Quantize ISF parameters (46 bits) every 40 ms                *
     *---------------------------------------------------------------*/
		/* Convert isps to frequency domain 0..6400 */
		E_LPC_isp_isf_conversion( &isp[( k + 2 ) * M], isfnew, M );
		/* quantize 1st and 2nd LPC with 46 bits */
		qpisf_2s_46b( isfnew, isfnew, past_isfq_2, prm_tcx, 4 );
		/* Convert isfs to the cosine domain */
		isf2isp( isfnew, ispnew_q, M );
		/* interpol quantized lpc */
		int_lpc_np1( &isp_q[k * M], ispnew_q, Aq, ( NB_SUBFR / 2 ), M );
		/*--------------------------------------------------*
     * Call medium TCX coder and find segmental SNR       *
     *--------------------------------------------------*/
		mvr2r( &synth[( k * L_DIV ) - M], synth_tcx, M );
		mem_w0_tcx = mem_w0[k];
		mem_xnq_tcx = mem_xnq[k];
		mvr2r( &wovlp[k * 128], wovlp_tcx, 128 );
		coder_tcx( Aq,
		    &speech[k * L_DIV],
		    &mem_w0_tcx,
		    &mem_xnq_tcx,
		    &synth_tcx[M],
		    exc_tcx,
		    wovlp_tcx,
		    ovlp_size[k],
		    2 * L_DIV,
		    ( ( nbits / 2 ) - NBITS_LPC ),
		    prm_tcx + NPRM_LPC );
		mem_wsyn_tcx = mem_wsyn[k];
		{
			float buf[L_FRAME_PLUS];
			find_wsp( &A[i1 * ( NB_SUBFR / 2 ) * ( M + 1 )], &synth_tcx[M], buf, &mem_wsyn_tcx, 2 * L_DIV );
			tmp = segsnr( &wsp[k * L_DIV], buf, 2 * L_DIV, L_SUBFR );
		}
		/*--------------------------------------------------------*
     * Save tcx parameters if tcx segmental SNR is better     *
     *--------------------------------------------------------*/
		if( tmp > snr1 ) {
			snr1 = tmp;
			for( i = 0; i < 2; i++ ) {
				mod[k + i] = 2;
				coding_mod[k + i] = 2;
			}
			mvr2r( ispnew_q, &isp_q[( k + 2 ) * M], M );
			mem_w0[k + 2] = mem_w0_tcx;
			mem_xnq[k + 2] = mem_xnq_tcx;
			mem_wsyn[k + 2] = mem_wsyn_tcx;
			ovlp_size[k + 2] = 64;
			mvr2r( past_isfq_2, past_isfq_1, M );
			mvr2r( wovlp_tcx, &wovlp[( k + 2 ) * 128], 128 );
			mvr2r( &synth_tcx[M], &synth[k * L_DIV], 2 * L_DIV );
			mvr2r( exc_tcx, &exc[k * L_DIV], 2 * L_DIV );
			mvi2i( prm_tcx, prm, NPRM_LPC + NPRM_TCX40 );
			/* lpc coefficient needed for HF extension */
			mvr2r( Aq, &AqLF[k * 4 * ( M + 1 )], 9 * ( M + 1 ) );
		}
		snr2 += 0.5f * snr1;
	} /* end of i1 */
	k = 0;
	/* set pointer to parameters */
	prm = param + ( k * NPRM_DIV );
	/*---------------------------------------------------------------*
   *  Quantize ISF parameters (46 bits) every 80 ms                *
   *---------------------------------------------------------------*/
	mvr2r( st->past_isfq, past_isfq_2, M );
	/* Convert isps to frequency domain 0..6400 */
	E_LPC_isp_isf_conversion( &isp[( k + 4 ) * M], isfnew, M );
	/* quantize 1st and 2nd LPC with 46 bits */
	qpisf_2s_46b( isfnew, isfnew, past_isfq_2, prm_tcx, 4 );
	/* Convert isfs to the cosine domain */
	isf2isp( isfnew, ispnew_q, M );
	/* interpol quantized lpc */
	int_lpc_np1( &isp_q[k * M], ispnew_q, Aq, NB_SUBFR, M );
	/*--------------------------------------------------*
   * Call long TCX coder and find segmental SNR       *
   *--------------------------------------------------*/
	mvr2r( &synth[( k * L_DIV ) - M], synth_tcx, M );
	mem_w0_tcx = mem_w0[k];
	mem_xnq_tcx = mem_xnq[k];
	mvr2r( &wovlp[k * 128], wovlp_tcx, 128 );
	coder_tcx( Aq,
	    &speech[k * L_DIV],
	    &mem_w0_tcx,
	    &mem_xnq_tcx,
	    &synth_tcx[M],
	    exc_tcx,
	    wovlp_tcx,
	    ovlp_size[k],
	    4 * L_DIV,
	    ( nbits - NBITS_LPC ),
	    prm_tcx + NPRM_LPC );
	mem_wsyn_tcx = mem_wsyn[k];
	{
		float buf[L_FRAME_PLUS];
		find_wsp( &A[0 * ( NB_SUBFR / 2 ) * ( M + 1 )], &synth_tcx[M], buf, &mem_wsyn_tcx, 4 * L_DIV );
		tmp = segsnr( &wsp[k * L_DIV], buf, 4 * L_DIV, L_SUBFR );
	}
	/*--------------------------------------------------------*
   * Save tcx parameters if tcx segmental SNR is better     *
   *--------------------------------------------------------*/
	if( tmp > snr2 ) {
		snr2 = tmp;
		for( i = 0; i < 4; i++ ) {
			mod[k + i] = 3;
			coding_mod[k + i] = 3;
		}
		mvr2r( ispnew_q, &isp_q[( k + 4 ) * M], M );
		mem_w0[k + 4] = mem_w0_tcx;
		mem_xnq[k + 4] = mem_xnq_tcx;
		mem_wsyn[k + 4] = mem_wsyn_tcx;
		ovlp_size[k + 4] = 128;
		mvr2r( past_isfq_2, past_isfq_1, M );
		mvr2r( wovlp_tcx, &wovlp[( k + 4 ) * 128], 128 );
		mvr2r( &synth_tcx[M], &synth[k * L_DIV], 4 * L_DIV );
		mvr2r( exc_tcx, &exc[k * L_DIV], 4 * L_DIV );
		mvi2i( prm_tcx, prm, NPRM_LPC + NPRM_TCX80 );
		/* lpc coefficient needed for HF extension */
		mvr2r( Aq, &AqLF[k * 4 * ( M + 1 )], 17 * ( M + 1 ) );
	}
	/*--------------------------------------------------*
   * Update filter memory.                            *
   *--------------------------------------------------*/
	st->old_ovlp_size = ovlp_size[NB_DIV];
	st->old_mem_w0 = mem_w0[NB_DIV];
	st->old_mem_xnq = mem_xnq[NB_DIV];
	st->old_mem_wsyn = mem_wsyn[NB_DIV];
	mvr2r( &wovlp[4 * 128], st->old_wovlp, 128 );
	mvr2r( past_isfq_1, st->past_isfq, M );
	mvr2r( &isp_q[NB_DIV * M], st->ispold_q, M );
	/*--------------------------------------------------*
   * Update exc for next frame.                       *
   *--------------------------------------------------*/
	mvr2r( &old_exc[L_FRAME_PLUS], st->old_exc, PIT_MAX_MAX + L_INTERPOL );
	mvr2r( old_d_wsp, st->old_d_wsp, PIT_MAX_MAX / OPL_DECIM ); /* d_wsp already shifted */
	return;
}
/*--------------------------------------------------*
 * Find weighted speech (formula from AMR-WB)       *
 *--------------------------------------------------*/
void find_wsp( float A[],
    float            speech[], /* speech[-M..lg]   */
    float            wsp[],    /* wsp[0..lg]       */
    float *          mem_wsp,  /* memory           */
    int              lg )
{
	int    i_subfr;
	float *p_A, Ap[M + 1];
	p_A = A;
	for( i_subfr = 0; i_subfr < lg; i_subfr += L_SUBFR ) {
		E_LPC_a_weight( p_A, Ap, GAMMA1, M );
		E_UTIL_residu( Ap, &speech[i_subfr], &wsp[i_subfr], L_SUBFR );
		p_A += ( M + 1 );
	}
	E_UTIL_deemph( wsp, TILT_FAC, lg, mem_wsp );
	return;
}
