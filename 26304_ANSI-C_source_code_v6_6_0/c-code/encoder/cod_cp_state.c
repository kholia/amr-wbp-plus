#include "../include/mem.h"
#include "../include/options.h"
#include "../include/proto_func.h"
#include "../lib_amr/enc_main.h"
#include "../lib_amr/typedef.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef WB_enc_if_state
typedef struct
{
	Word16 sid_update_counter; /* Number of frames since last SID */
	Word16 sid_handover_debt;  /* Number of extra SID_UPD frames to schedule */
	Word16 prev_ft;            /* Type of the previous frame */
	void * encoder_state;      /* Points encoder state structure */
} WB_enc_if_state;
#endif
#define Q_MAX 8 /* scaling max for signal                 */
void E_MAIN_reset( void *st, Word16 reset_all );
void copy_coder_state(
    Coder_State_Plus *wbP, /* AMR-WB+ state struct */
    void *            st,  /* AMR-WB state struct  */
    short             sw,  /* sw=0 -> switch from WB to WB+, sw=1 -> switch from WB+ to WB */
    short             use_case_mode )
{
	Word16           Q_new, tmp, tmp2;
	Word32           i;
	float            ftemp;
	WB_enc_if_state *s1;
	Coder_State *    wb;
	s1 = (WB_enc_if_state *)st;
	wb = (Coder_State *)( s1->encoder_state );
	if( sw == 0 ) {
		short t;
		init_coder_amrwb_plus( wbP, 1, 0, use_case_mode, 1 );
		mvr2r( wb->mem_isf, wbP->isfold, M );
		mvs2r( wb->mem_isf_q, wbP->past_isfq, M );
		mvr2r( wb->mem_isp, wbP->ispold, M );
		for( i = 0; i < M; i++ ) {
			wbP->ispold_q[i] = (float)( wb->mem_isp_q[i] * 1024 );
			wbP->ispold_q[i] = wbP->ispold_q[i] / ( 1024 * 32767 );
		}
		mvi2i( (int *)wb->mem_ol_lag, (int *)wbP->old_ol_lag, 5 );
		wbP->old_T0_med = (int)( wb->mem_T0_med );
		mvr2r( wb->mem_decim, wbP->right.mem_decim, 30 );
		mvr2r( wb->mem_sig_in, wbP->right.mem_sig_in, 4 );
		wbP->right.mem_preemph = wb->mem_preemph;
		mvr2r( wb->mem_speech + 48, wbP->old_speech_pe, 16 + 512 );
		t = PIT_MAX_MAX; //- PIT_MAX;
		mvr2r( wb->mem_wsp, wbP->old_d_wsp + ( PIT_MAX_MAX / OPL_DECIM ) - ( PIT_MAX / OPL_DECIM ), PIT_MAX / OPL_DECIM );
		mvr2r( wb->mem_hp_wsp, wbP->hp_old_wsp + ( PIT_MAX_MAX / OPL_DECIM ) - ( PIT_MAX / OPL_DECIM ), PIT_MAX / OPL_DECIM );
		wbP->mem_wsp = wb->mem_wsp_df;
		mvr2r( wb->mem_decim2, wbP->mem_lp_decim2, 3 );
		wbP->old_T0_med = wb->mem_T0_med;
		wbP->ol_gain = wb->mem_ol_gain;
		mvr2r( wb->mem_hf_wsp + 2, wbP->hp_ol_ltp_mem, 7 );
		wbP->ada_w = wb->mem_ada_w;
		wbP->ol_wght_flg = (short)wb->mem_ol_wght_flg;
		wbP->old_mem_w0 = wb->mem_w0;
		mvr2r( wb->mem_syn, wbP->old_synth, M );
		for( i = 0; i < PIT_MAX + L_INTERPOL; i++ ) {
			wbP->old_exc[PIT_MAX_MAX - PIT_MAX + i] = (float)( wb->mem_exc[i] * pow( 2, -( wb->mem_q ) ) );
		}
		wbP->old_ovlp_size = 0;
		wbP->SwitchFlagPlusToWB = 1;
		if( use_case_mode == USE_CASE_B ) {
			mvr2r( wb->vadSt->mem_bckr_est, wbP->vadSt->bckr_est, 12 );
			mvr2r( wb->vadSt->mem_ave_level, wbP->vadSt->ave_level, 12 );
			mvr2r( wb->vadSt->mem_level, wbP->vadSt->old_level, 12 );
			mvr2r( wb->vadSt->mem_sub_level, wbP->vadSt->sub_level, 12 );
			mvr2r( &( wb->vadSt->mem_a_data5[0][0] ), &( wbP->vadSt->a_data5[0][0] ), 10 );
			mvr2r( wb->vadSt->mem_a_data3, wbP->vadSt->a_data3, 6 );
			wbP->vadSt->sp_max = wb->vadSt->mem_sp_max;
			wbP->vadSt->speech_level = wb->vadSt->mem_speech_level;
			wbP->vadSt->burst_count = wb->vadSt->mem_burst_count;
			wbP->vadSt->hang_count = wb->vadSt->mem_hang_count;
			wbP->vadSt->stat_count = wb->vadSt->mem_stat_count;
			wbP->vadSt->vadreg = wb->vadSt->mem_vadreg;
			wbP->vadSt->pitch_tone = wb->vadSt->mem_pitch_tone;
			wbP->vadSt->sp_est_cnt = wb->vadSt->mem_sp_est_cnt;
			wbP->vadSt->sp_max_cnt = wb->vadSt->mem_sp_max_cnt;
			if( wbP->stClass->StatClassCount == 0 ) {
				wbP->stClass->StatClassCount = 15;
			}
		}
	}
	else if( sw == 1 ) {
		E_MAIN_reset( (void *)wb, 1 );
		mvr2r( wbP->isfold, wb->mem_isf, M );
		mvr2s( wbP->past_isfq, wb->mem_isf_q, M );
		mvr2r( wbP->ispold, wb->mem_isp, M );
		for( i = 0; i < M; i++ ) {
			wb->mem_isp_q[i] = (short)( wbP->ispold_q[i] * 32767 );
		}
		mvi2i( (int *)wbP->old_ol_lag, (int *)wb->mem_ol_lag, 5 );
		wb->mem_T0_med = ( Word32 )( wbP->old_T0_med );

		mvr2r( wbP->right.mem_decim, wb->mem_decim, 30 );
		mvr2r( wbP->right.mem_sig_in, wb->mem_sig_in, 4 );
		wb->mem_preemph = wbP->right.mem_preemph;
		mvr2r( wbP->old_speech_pe, wb->mem_speech + 48, 16 + 512 );

		for( i = 0; i < 48; i++ )
			wb->mem_speech[i] = 0.0f; // done in the reset

		mvr2r( wbP->old_d_wsp + ( PIT_MAX_MAX / OPL_DECIM ) - ( PIT_MAX / OPL_DECIM ), wb->mem_wsp, PIT_MAX / OPL_DECIM );
		mvr2r( wbP->hp_old_wsp + ( PIT_MAX_MAX / OPL_DECIM ) - ( PIT_MAX / OPL_DECIM ), wb->mem_hp_wsp, PIT_MAX / OPL_DECIM );
		wb->mem_wsp_df = wbP->mem_wsp;
		mvr2r( wbP->mem_lp_decim2, wb->mem_decim2, 3 );
		wb->mem_ol_gain = wbP->ol_gain;
		mvr2r( wbP->hp_ol_ltp_mem, wb->mem_hf_wsp + 2, 7 );
		wb->mem_hf_wsp[0] = 0.0f;
		wb->mem_hf_wsp[1] = 0.0f;
		wb->mem_ada_w = wbP->ada_w;
		wb->mem_ol_wght_flg = (UWord8)wbP->ol_wght_flg;
		wb->mem_w0 = wbP->old_mem_w0;
		if( wbP->prev_mod == 0 ) {
			for( i = 0; i < 4; i++ ) {
				ftemp = 20.0f * (float)log10( (float)wbP->mem_gain_code[3 - i] ) * 1024.0f;
				if( ftemp < -32768.0f )
					wb->mem_gain_q[i] = -32768;
				else if( ftemp > 32767.0f )
					wb->mem_gain_q[i] = 32767;
				else
					wb->mem_gain_q[i] = (Word16)ftemp;
			}
		}
		mvr2r( wbP->old_synth, wb->mem_syn, M );
		mvr2r( wbP->old_synth, wb->mem_syn2, M );
		Q_new = 0;
		wb->mem_q = 0;
		tmp = ( Word16 )( fabs( wbP->old_exc[0] ) );
		for( i = 0; i < ( PIT_MAX + L_INTERPOL ); i++ ) {
			tmp2 = ( Word16 )( fabs( wbP->old_exc[PIT_MAX_MAX - PIT_MAX + i] ) );
			if( tmp < tmp2 )
				tmp = tmp2;
		}
		while( ( tmp < 0x4000 ) && ( Q_new < 3 ) ) {
			tmp = ( tmp << 1 );
			Q_new = Q_new + 1;
		}
		for( i = 0; i < ( PIT_MAX + L_INTERPOL ); i++ ) {
			wb->mem_exc[i] = ( ( Word16 )( wbP->old_exc[PIT_MAX_MAX - PIT_MAX + i] ) ) << Q_new;
		}
		wb->mem_q = (Word16)Q_new;
		for( i = 0; i < 4; i++ )
			wb->mem_subfr_q[i] = (Word16)Q_new;
		if( use_case_mode == USE_CASE_B ) {
			mvr2r( wbP->vadSt->bckr_est, wb->vadSt->mem_bckr_est, 12 );
			mvr2r( wbP->vadSt->ave_level, wb->vadSt->mem_ave_level, 12 );
			mvr2r( wbP->vadSt->old_level, wb->vadSt->mem_level, 12 );
			mvr2r( wbP->vadSt->sub_level, wb->vadSt->mem_sub_level, 12 );
			mvr2r( &( wbP->vadSt->a_data5[0][0] ), &( wb->vadSt->mem_a_data5[0][0] ), 10 );
			mvr2r( wbP->vadSt->a_data3, wb->vadSt->mem_a_data3, 6 );
			wb->vadSt->mem_sp_max = wbP->vadSt->sp_max;
			wb->vadSt->mem_speech_level = wbP->vadSt->speech_level;
			wb->vadSt->mem_burst_count = wbP->vadSt->burst_count;
			wb->vadSt->mem_hang_count = wbP->vadSt->hang_count;
			wb->vadSt->mem_stat_count = wbP->vadSt->stat_count;
			wb->vadSt->mem_vadreg = wbP->vadSt->vadreg;
			wb->vadSt->mem_pitch_tone = wbP->vadSt->pitch_tone;
			wb->vadSt->mem_sp_est_cnt = wbP->vadSt->sp_est_cnt;
			wb->vadSt->mem_sp_max_cnt = wbP->vadSt->sp_max_cnt;
		}
	}
	else {
		fprintf( stderr, " Unknown mode switching paramer: %d \n", sw );
		exit( -1 );
	}
	return;
}
