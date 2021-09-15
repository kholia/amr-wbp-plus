#ifndef _stereo_cb_lo_
#define _stereo_cb_lo_

#define LO_NUMSTAGES_CB_LO_1500 3
#define LO_NUMSTAGES_CB_LO_1900 4
#define LO_NUMSTAGES_CB_LO_2500 5
#define LO_NUMSTAGES_CB_LO_3000 6
// low band MA-SMSVQ means for MA prediction 16 and 24 kHz
extern const float stereo_lo_mean[VDIM];
extern const float stereo_lo_smsvq_a_CB_LO_1500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_b_CB_LO_1500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_c_CB_LO_1500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_a_CB_LO_1900[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_b_CB_LO_1900[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_c_CB_LO_1900[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_d_CB_LO_1900[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_a_CB_LO_2500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_b_CB_LO_2500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_c_CB_LO_2500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_d_CB_LO_2500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_e_CB_LO_2500[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_a_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_b_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_c_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_d_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_e_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_lo_smsvq_f_CB_LO_3000[LO_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
#endif
