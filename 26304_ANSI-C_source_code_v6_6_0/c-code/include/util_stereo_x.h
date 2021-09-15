#ifndef _UTIL_STEREO_X_H
#define _UTIL_STEREO_X_H

typedef struct
{
	const int     vdim;
	const int     nstages;
	const int     intens;
	const int *   cbsizes;
	const float **cbs;
} MSVQ;

typedef struct
{
	const float  a;    /* Prediction coefficient */
	const float  a_fe; /* Prediction coefficient to be applied in decoder in case of fe */
	const float *mean; /* Mean vector */
	const MSVQ   msvq;
} PMSVQ;

#endif
