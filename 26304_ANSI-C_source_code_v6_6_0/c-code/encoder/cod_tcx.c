#include "../include/amr_plus.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void coder_tcx(
    float  A[],       /* input:  coefficients NxAz[M+1]  */
    float  speech[],  /* input:  speech[-M..lg]          */
    float *mem_wsp,   /* in/out: memory of target        */
    float *mem_wsyn,  /* in/out: memory of quantized xn  */
    float  synth[],   /* in/out: synth[-M..lg]           */
    float  exc[],     /* output: exc[0..lg]              */
    float  wovlp[],   /* i/o:    wovlp[0..127]           */
    int    ovlp_size, /* input:  0, 64 or 128 (0=acelp)  */
    int    L_frame,   /* input:  frame length            */
    int    nb_bits,   /* input:  number of bits allowed  */
    int    prm[] )       /* output: tcx parameters          */
{
	int    i, i_subfr, lext, lg, index;
	float  tmp, gain, fac_ns;
	float *p_A, Ap[M + 1];
	float *xri, *xn, *xnq;
	float  window[256];
	xn = synth;
	xnq = xri = exc;
	/*------ set length of overlap (lext) and length of encoded frame (lg) -----*/
	lext = L_OVLP;
	if( L_frame == ( L_FRAME_PLUS / 2 ) ) {
		lext = L_OVLP / 2;
	}
	if( L_frame == ( L_FRAME_PLUS / 4 ) ) {
		lext = L_OVLP / 4;
	}
	lg = L_frame + lext;
	/* built window for overlaps section */
	cos_window( window, ovlp_size, lext );
	/*-----------------------------------------------------------*
   * Find target xn[] (weighted speech when prev mode is TCX)  *
   * Note that memory isn't updated in the overlap area.       *
   *-----------------------------------------------------------*/
	p_A = A;
	for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR ) {
		E_LPC_a_weight( p_A, Ap, GAMMA1, M );
		E_UTIL_residu( Ap, &speech[i_subfr], &xn[i_subfr], L_SUBFR );
		p_A += ( M + 1 );
	}
	E_UTIL_deemph( xn, TILT_FAC, L_frame, mem_wsp );
	/* overlap area (xn[L_frame..L_frame+lext]) */
	E_LPC_a_weight( p_A, Ap, GAMMA1, M );
	E_UTIL_residu( Ap, &speech[L_frame], &xn[L_frame], lext );
	tmp = *mem_wsp;
	E_UTIL_deemph( &xn[L_frame], TILT_FAC, lext, &tmp );
	/* remove weighted ZIR when previous frame is ACELP */
	if( ovlp_size == 0 ) {
		for( i = 0; i < ( 2 * L_SUBFR ); i++ ) {
			xn[i] -= wovlp[i];
		}
	}
	/* xn[] windowing for TCX overlap */
	for( i = 0; i < ovlp_size; i++ ) {
		xn[i] *= window[i];
	}
	for( i = 0; i < lext; i++ ) {
		xn[L_frame + i] *= window[ovlp_size + i];
	}
	/*-----------------------------------------------------------*
  * Compute the FFT of xn[].                                  *
  * Coefficients (xri[]) order are                            *
  *    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
  * Note that last FFT element (re[n/2]) is zeroed.           *
  *-----------------------------------------------------------*/
	fft9( xn, xri, (short)lg );
	xri[1] = 0.0; /* freq bin at 6400 Hz zeroed */
	              /*-----------------------------------------------------------*
  * Spectral algebraic quantization                           *
  * with adaptive low frequency emphasis/deemphasis.          *
  * Noise factor is the average level of unquantized freq.    *
  *-----------------------------------------------------------*/
	adap_low_freq_emph( xri, lg );
	nb_bits -= ( 3 + 7 ); /* fac_ns = 3 bits, gain = 7 bits */
	/* remove also the redundancy bits of the TCX gain
     TCX-40 -> 6 bits (in 2nd packet)
     TCX-80 -> 9 bits (3 bits in 2nd, 3rd and 4th packet) */
	if( L_frame == ( L_FRAME_PLUS / 2 ) ) {
		nb_bits -= 6;
	}
	if( L_frame == L_FRAME_PLUS ) {
		nb_bits -= 9;
	}
	fac_ns = AVQ_cod( xri, prm + 2, nb_bits, lg / 8 );
	for( i = 0; i < lg; i++ ) {
		xri[i] = (float)prm[i + 2];
	}
	adap_low_freq_deemph( xri, lg );
	/* quantize noise factor (noise factor = 0.1 to 0.8) */
	tmp = 8.0f - ( 10.0f * fac_ns );
	index = (int)floor( tmp + 0.5 );
	if( index < 0 ) {
		index = 0;
	}
	if( index > 7 ) {
		index = 7;
	}
	prm[0] = index; /* fac_ns : 3 bits */
	                /*-----------------------------------------------------------*
  * Compute inverse FFT for obtaining xnq[] without noise.    *
  * Coefficients (xri[]) order are                            *
  *    re[0], re[n/2], re[1], re[2], ... re[n/2-1], im[n/2-1] *
  * Note that last FFT element (re[n/2]) is zeroed.           *
  *-----------------------------------------------------------*/
	xri[1] = 0.0;   /* freq bin at 6400 Hz zeroed */
	ifft9( xri, xnq, (short)lg );
	/*-----------------------------------------------------------*
  * find and quantize gain, multiply xnq[] by gain.           *
  * windowing of xnq[] for TCX overlap.                       *
  *-----------------------------------------------------------*/
	gain = get_gain( xn, xnq, lg );
	prm[1] = q_gain_tcx( xnq, lg, &gain );
	for( i = 0; i < lg; i++ ) {
		xnq[i] *= gain;
	}
	/* adaptive windowing on overlap (beginning and end of frame) */
	for( i = 0; i < ovlp_size; i++ ) {
		xnq[i] *= window[i];
	}
	for( i = 0; i < lext; i++ ) {
		xnq[i + L_frame] *= window[ovlp_size + i];
	}
	/*-----------------------------------------------------------*
  * TCX overlap and add.  Update memory for next overlap.     *
  *-----------------------------------------------------------*/
	for( i = 0; i < L_OVLP; i++ ) {
		xnq[i] += wovlp[i];
	}
	/* save overlap for next frame */
	for( i = 0; i < lext; i++ ) {
		wovlp[i] = xnq[i + L_frame];
	}
	for( i = lext; i < L_OVLP; i++ ) {
		wovlp[i] = 0.0;
	}
	/*-----------------------------------------------------------*
  * find excitation and synthesis                             *
  *-----------------------------------------------------------*/
	E_UTIL_f_preemph( xnq, TILT_FAC, L_frame, mem_wsyn );
	p_A = A;
	for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR ) {
		E_LPC_a_weight( p_A, Ap, GAMMA1, M );
		E_UTIL_synthesis( Ap, &xnq[i_subfr], &synth[i_subfr], L_SUBFR, &synth[i_subfr - M], 0 );
		E_UTIL_residu( p_A, &synth[i_subfr], &exc[i_subfr], L_SUBFR );
		p_A += ( M + 1 );
	}
	return;
}
