#ifndef _COD_HI_STEREO
#define _COD_HI_STEREO

void cod_hi_stereo( float speech_hi[],
    float                 right_hi[],
    float                 AqLF[],
    int                   prm[],
    Coder_State_Plus *    st );

void dec_hi_stereo( float synth_hi_t0[],
    float                 right_hi[],
    float                 left_hi[],
    float                 AqLF[],
    int                   prm[],
    int                   bad_frame[],
    int                   fscale,
    Decoder_State_Plus *  st );

void init_cod_hi_stereo( Coder_State_Plus *st );
void init_dec_hi_stereo( Decoder_State_Plus *st );
#endif
