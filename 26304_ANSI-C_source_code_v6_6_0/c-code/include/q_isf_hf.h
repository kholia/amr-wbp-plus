#ifndef q_isf_hf_h
#define q_isf_hf_h
/*
*/

/*-------------------------------------------------------------------*
 *                         q_isf_hf.h
 *-------------------------------------------------------------------*
 * Quantization tables for two-stage of HF ISFs
 *-------------------------------------------------------------------*/

/*
Use speech (57 min) and audio (31 min) database

result from kmean v05:

mu = 0.500000,	scale = 1.000000,	model = AR
frames counter = 206830

average spectral distortion =  1.433384
(0.0-0.5] dB spectral distortion =  0.287676
(0.5-1.0] dB spectral distortion =  21.727506
(1.0-1.5] dB spectral distortion =  41.181163
(1.5-2.0] dB spectral distortion =  22.323647
(2.0-2.5] dB spectral distortion =  9.506358
(2.5-3.0] dB spectral distortion =  3.369434
(3.0-3.5] dB spectral distortion =  1.058841
(3.5-4.0] dB spectral distortion =  0.354881
(4.0-...] dB spectral distortion =  0.190495
*/

#define Q_ISF_ORDER 8 /* order of linear prediction filter */

#define SIZE_BK1_HF 4
#define SIZE_BK2_HF 128

/* Actual data in tables_plus.c */

#endif
