#include "../include/mem.h"
#include "../include/proto_func.h"
#include "../lib_amr/dec_main.h"
#include "../lib_amr/typedef.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef WB_dec_if_state
typedef struct
{
	Word16 reset_flag_old; /* previous was homing frame  */
	Word16 prev_ft;        /* previous frame type        */
	Word16 prev_mode;      /* previous mode              */
	void * decoder_state;  /* Points decoder state       */
} WB_dec_if_state;
#endif
#define Q_MAX 8 /* scaling max for signal                 */
void D_MAIN_reset( void *st, Word16 reset_all );
void copy_decoder_state(
    Decoder_State_Plus *wbP,
    void *              st,
    short               sw )
{
	Word16           Q_new, tmp, tmp2;
	Word32           i, ltemp;
	float            ftemp;
	WB_dec_if_state *s1;
	Decoder_State *  wb;
	s1 = (WB_dec_if_state *)st;
	wb = (Decoder_State *)( s1->decoder_state );
	if( sw == 0 ) {
		init_decoder_amrwb_plus( wbP, 1, 0, 1 );
		for( i = 0; i < M; i++ ) {
			wbP->ispold[i] = (float)wb->mem_isp[i] / 32768.0f;
			wbP->isfold[i] = (float)( wb->mem_isf[i] );
			wbP->past_isfq[i] = (float)( wb->mem_isf_q[i] );
		}
		for( i = 0; i < L_MEANBUF; i++ ) {
			mvr2r( (float *)wbP->isfold, &( wbP->isf_buf[i * M] ), M );
		}
		wbP->old_T0 = (int)( wb->mem_T0 );
		wbP->old_T0_frac = (int)( wb->mem_T0_frac );
		wbP->old_T_pf[0] = wbP->old_T_pf[1] = wbP->old_T0;
		wbP->old_gain_pf[0] = wbP->old_gain_pf[1] = 0.0f;
		for( i = 0; i < M; i++ ) {
			ltemp = ( wb->mem_syn_hi[i] << 12 ) + wb->mem_syn_lo[i];
			ltemp = ltemp >> 8;
			wbP->old_synth[i] = (float)ltemp;
		}
		for( i = 0; i < PIT_MAX + L_INTERPOL; i++ ) {
			wbP->old_exc[i + PIT_MAX_MAX - PIT_MAX] = (float)( wb->mem_exc[i] * pow( 2, -( wb->mem_q ) ) );
		}
		wbP->mem_deemph = (float)wb->mem_deemph;
		wbP->mem_sig_out[0] = (float)wb->mem_sig_out[2] * 2.0f;
		wbP->mem_sig_out[1] = (float)wb->mem_sig_out[0] * 2.0f;
		wbP->mem_sig_out[2] = (float)wb->mem_sig_out[4];
		wbP->mem_sig_out[3] = (float)wb->mem_sig_out[5];
		for( i = 0; i < 2 * L_FILT; i++ ) {
			wbP->right.mem_oversamp[i] = (float)wb->mem_oversamp[i];
		}
		for( i = 0; i < PIT_MAX + L_SUBFR; i++ ) {
			wbP->old_synth_pf[PIT_MAX_MAX - PIT_MAX + i] = (float)wb->mem_syn_out[i];
		}
		wbP->gc_threshold = (float)wb->mem_gc_thres * (float)pow( 2.0f, -16.0f );
		wbP->ramp_state = 0;
	}
	else {
		D_MAIN_reset( (void *)wb, 1 );
		for( i = 0; i < M; i++ ) {
			wb->mem_isp[i] = ( Word16 )( wbP->ispold[i] * 32768.0f );
			wb->mem_isf[i] = ( Word16 )( wbP->isfold[i] );
			wb->mem_isf_q[i] = ( Word16 )( wbP->past_isfq[i] );
		}
		wb->mem_T0 = ( Word16 )( wbP->old_T0 );
		wb->mem_T0_frac = ( Word16 )( wbP->old_T0_frac );
		for( i = 0; i < M; i++ ) {
			ltemp = ( Word32 )( wbP->old_synth[i] * 256.0f );
			wb->mem_syn_hi[i] = ( Word16 )( ltemp >> 12 );
			wb->mem_syn_lo[i] = ( Word16 )( ltemp & 0x00000FFF );
		}
		if( wbP->last_mode == 0 ) {
			for( i = 0; i < 4; i++ ) {
				ftemp = 20.0f * (float)log10( (float)wbP->mem_gain_code[3 - i] ) * 1024.0f;
				if( ftemp < -32768.0f )
					wb->mem_gain[i] = (Word16)-32768;
				else if( ftemp > 32767.0f )
					wb->mem_gain[i] = (Word16)32767;
				else
					wb->mem_gain[i] = (Word16)ftemp;
			}
		}
		Q_new = 0;
		wb->mem_q = 0;
		tmp = ( Word16 )( fabs( wbP->old_exc[0] ) );
		for( i = 0; i < ( PIT_MAX + L_INTERPOL ); i++ ) {
			tmp2 = ( Word16 )( fabs( wbP->old_exc[i + PIT_MAX_MAX - PIT_MAX] ) );
			if( tmp < tmp2 )
				tmp = tmp2;
		}
		while( ( tmp < 0x4000 ) && ( Q_new < 3 ) ) {
			tmp = ( tmp << 1 );
			Q_new = Q_new + 1;
		}
		for( i = 0; i < ( PIT_MAX + L_INTERPOL ); i++ ) {
			wb->mem_exc[i] = ( ( Word16 )( wbP->old_exc[i + PIT_MAX_MAX - PIT_MAX] ) ) << Q_new;
		}
		wb->mem_q = (Word16)Q_new;
		for( i = 0; i < 4; i++ ) {
			wb->mem_subfr_q[i] = (Word16)Q_new;
		}
		wb->mem_deemph = (Word16)wbP->mem_deemph;
		wb->mem_sig_out[0] = ( Word16 )( wbP->mem_sig_out[1] * 0.5f );
		wb->mem_sig_out[1] = (Word16)0;
		wb->mem_sig_out[2] = ( Word16 )( wbP->mem_sig_out[0] * 0.5f );
		wb->mem_sig_out[3] = (Word16)0;
		wb->mem_sig_out[4] = ( Word16 )( wbP->mem_sig_out[2] );
		wb->mem_sig_out[5] = ( Word16 )( wbP->mem_sig_out[3] );
		for( i = 0; i < 2 * L_FILT; i++ ) {
			wb->mem_oversamp[i] = (Word16)wbP->right.mem_oversamp[i];
		}
		for( i = 0; i < PIT_MAX + L_SUBFR; i++ ) {
			ftemp = wbP->old_synth_pf[PIT_MAX_MAX - PIT_MAX + i];
			if( ftemp < -32768.0f )
				wb->mem_syn_out[i] = (Word16)-32768;
			else if( ftemp > 32767.0f )
				wb->mem_syn_out[i] = (Word16)32767;
			else
				wb->mem_syn_out[i] = (Word16)ftemp;
		}
		wb->mem_gc_thres = ( Word32 )( wbP->gc_threshold * (float)pow( 2.0f, 16.0f ) + 0.5f );
		for( i = 0; i < 9; i++ )
			wb->lpc_hf_plus[i] = wbP->mem_lpc_hf[i];
		wb->gain_hf_plus = wbP->mem_gain_hf;
		for( i = 0; i < 2 * L_FILT; i++ )
			wb->mem_oversamp_hf_plus[i] = wbP->right.mem_oversamp_hf[i];
		for( i = 0; i < 8; i++ )
			wb->mem_syn_hf_plus[i] = wbP->right.mem_syn_hf[i];
		wb->threshold_hf = wbP->right.threshold;
		wb->lp_amp_hf = wbP->right.lp_amp;
		wb->ramp_state = 64;
	}
	return;
}
