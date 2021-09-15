#ifndef _TCX_STEREO_H
#define _TCX_STEREO_H

void init_tcx_stereo_encoder( Coder_State_Plus *st );
void init_tcx_stereo_decoder( Decoder_State_Plus *st );

void cod_tcx_stereo( float mono_2k[],
    float                  right_2k[],
    int                    param[],
    int                    brMode,
    int                    mod[],
    int                    fscale,
    Coder_State_Plus *     st );

void dec_tcx_stereo( float synth_2k[],
    float                  left_2k[],
    float                  right_2k[],
    int                    param[],
    int                    bad_frame[],
    Decoder_State_Plus *   st );

#endif
