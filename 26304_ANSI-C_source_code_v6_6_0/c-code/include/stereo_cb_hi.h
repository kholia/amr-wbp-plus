
#ifndef _stereo_cb_hi_
#define _stereo_cb_hi_

//#define HI_NUMSTAGES_CB_HI_500 1
// Pat
#define HI_NUMSTAGES_CB_HI_500 2
#define HI_NUMSTAGES_CB_HI_1000 2
#define HI_NUMSTAGES_CB_HI_1500 3

extern const float stereo_hi_mean[VDIM];

extern const float stereo_hi_smsvq_a_CB_HI_500[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];

// Pat
extern const float stereo_hi_smsvq_b_CB_HI_500[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];

extern const float stereo_hi_smsvq_a_CB_HI_1000[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_hi_smsvq_b_CB_HI_1000[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];

extern const float stereo_hi_smsvq_a_CB_HI_1500[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_hi_smsvq_b_CB_HI_1500[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
extern const float stereo_hi_smsvq_c_CB_HI_1500[HI_SMSVQ_CBSIZE][( VDIM + 1 ) / 2];
#endif
