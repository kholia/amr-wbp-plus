#include "../include/amr_plus.h"
#include <math.h>
#ifndef PI2
#define PI2 6.283185307F
#endif
#define MAX_VECT_DIM 16
#define EPSILON 1e-16f
#define BUF_SIZ ( L_SUBFR + M )
float my_max( float x, float y )
{
	if( x > y ) {
		return x;
	}
	return y;
}
float my_min( float x, float y )
{
	if( x < y ) {
		return x;
	}
	return y;
}
void fir_filt( float *h, /* input : filter coefficients              */
    int               m, /* input : order of filter                  */
    float *           x, /* input : input signal (usually mono)      */
    float *           y, /* output: output signal (usually stereo)   */
    int               l  /* input : size of filtering                */
)
{
	float s;
	int   i, j;
	for( i = 0; i < l; i++ ) {
		s = 0.0f;
		for( j = 0; j < m; j++ ) {
			s += h[j] * x[i - j];
		}
		y[i] = s;
	}
	return;
}
void band_join_2k( float sig[],
    float                sig_2k[],
    float                sig_hi[],
    int                  lg )
{
	int i;
	/* interpolate to 12.8k fs*/
	interpol( sig_2k - L_FDEL_2k, sig, lg, (float *)filter_2k, L_FDEL_2k, 32, 5, 32.0f );
	for( i = 0; i < lg; i++ ) {
		sig[i] += sig_hi[i - L_FDEL];
	}
	return;
}
void band_split_taligned_2k( float sig[],
    float                          sig_2k[],
    float                          sig_hi[],
    int                            lg )
{
	int i;
	/* decimate to 2k fs */
	interpol( sig - L_FDEL, sig_2k, lg * 5 / 32, (float *)filter_2k, L_FDEL, 5, 32, 5.0f ); /*64 coeff*/
	/* interpolate to 12.8k fs*/
	interpol( sig_2k - L_FDEL_2k, sig_hi, lg, (float *)filter_2k, L_FDEL_2k, 32, 5, 32.0f ); /*10 coeff*/
	/* buf equal low frequency sig  + delay = 64 + 10*32/5 = 128 samples @ 12.8 kHz = 10ms delay*/
	for( i = 0; i < lg; i++ ) {
		sig_hi[i] = sig[i - 2 * L_FDEL] - sig_hi[i];
	}
	return;
}
void mix_ch( float *ch_left,  /* input: samples from left channel */
    float *         ch_right, /* input: samples from right channel */
    float *         ch_sum,   /* output: mixed mono signal */
    int             n,        /* input: length of frame */
    float           gA,       /* input: weight factor for left channel */
    float           gB        /* input: weight factor for right channel */
)
{
	int i;
	for( i = 0; i < n; i++ ) {
		ch_sum[i] = 0.5f * ( gA * ch_left[i] + gB * ch_right[i] );
	}
}

/************************************/
static int min_msvq( int a, /* (i)	:	Vector 1	*/
    int                  b  /* (i)	:	Vector 2	*/
)
{
	if( a < b ) {
		return a;
	}
	else {
		return b;
	}
}
/****************************************************************************/
static void m_cbcod( int mpfade[INTENS_MAX][MAX_NUMSTAGES], /* TBD	*/
    float *              dist,                              /* TBD	*/
    float *              x,                                 /* TBD	*/
    int *                pfad,                              /* TBD*/
    int                  stage,                             /* TBD*/
    const float *        cb,                                /* TBD*/
    int                  cbsize,                            /* TBD*/
    int                  m,                                 /* TBD*/
    int *                m_best,                            /* TBD*/
    float *              m_dist,                            /* TBD*/
    int                  vdim                               /*  -> vdim		: vector dimension */
)
{
	int          i, j;
	const float *y;
	float        distanz, d;
	for( i = 0; i < cbsize; i++ ) {
		y = &cb[i * vdim];
		distanz = 0;
		for( j = 0; j < vdim; j++ ) {
			d = x[j] - y[j];
			distanz += d * d;
			if( distanz >= *m_dist )
				goto endloop;
		}
		*m_dist = dist[*m_best] = distanz;
		for( j = 0; j < stage; j++ ) {
			mpfade[*m_best][j] = pfad[j];
		}
		mpfade[*m_best][stage] = i;
		for( j = 0; j < m; j++ ) {
			if( dist[j] > *m_dist ) {
				*m_dist = dist[j];
				*m_best = j;
			}
		}
	endloop:;
	}
}
/****************************************************************************/
static void msvq( float *y, int *inds, float *x, const MSVQ *msvq_tab )
/*  <- y			: reconstruction vector */
/*  <- inds		: the best path through the msvq trellis */
/*  -> x			: vector to be encoded */
/*  -> msvq_tab   : all required parameters and tables */
{
	int   i, j, k, l, pfadanz_max, arg;
	int   m_best;
	float m_dist;
	float e[MAX_VECT_DIM];
	int   mpfade_mem1[INTENS_MAX][MAX_NUMSTAGES];
	int   mpfade_mem2[INTENS_MAX][MAX_NUMSTAGES];
	int( *mpfade )[INTENS_MAX][MAX_NUMSTAGES];
	int( *alt_mpfade )[INTENS_MAX][MAX_NUMSTAGES];
	int( *tmp_pfad )[INTENS_MAX][MAX_NUMSTAGES];
	float         dist[MAX_NUMSTAGES][INTENS_MAX];
	const int *   cbsize = msvq_tab->cbsizes;
	const int     vdim = msvq_tab->vdim;
	const int     stages = msvq_tab->nstages;
	const int     m = msvq_tab->intens;
	const float **cb = msvq_tab->cbs;
	mpfade = &mpfade_mem1;
	alt_mpfade = &mpfade_mem2;
	for( j = 0; j < m; j++ ) {
		dist[0][j] = 1e30f;
	}
	m_best = 0;
	m_dist = 1e30f;
	pfadanz_max = cbsize[0];
	m_cbcod( *mpfade, dist[0], x, NULL, 0, cb[0], cbsize[0], m, &m_best, &m_dist, vdim );
	for( i = 1; i < stages; i++ ) {
		for( j = 0; j < m; j++ ) {
			dist[i][j] = 1e30f;
		}
		m_best = 0;
		m_dist = 1e30f;
		tmp_pfad = mpfade;
		mpfade = alt_mpfade;
		alt_mpfade = tmp_pfad;
		if( m < pfadanz_max ) {
			arg = m;
		}
		else {
			arg = pfadanz_max;
		}
		for( k = 0; k < min_msvq( m, pfadanz_max ); k++ ) {
			for( j = 0; j < vdim; j++ ) {
				e[j] = x[j];
				for( l = 0; l < i; l++ ) {
					e[j] -= cb[l][( *alt_mpfade )[k][l] * vdim + j];
				}
			}
			m_cbcod( *mpfade, dist[i], e, ( *alt_mpfade )[k], i, cb[i], cbsize[i], ( i == stages - 1 ? 1 : m ), &m_best, &m_dist, vdim );
		}
		pfadanz_max *= cbsize[i];
	}
	for( j = 0; j < stages; j++ ) {
		inds[j] = ( *mpfade )[0][j];
	}
	for( j = 0; j < vdim; j++ ) {
		y[j] = cb[0][inds[0] * vdim + j];
		for( l = 1; l < stages; l++ ) {
			y[j] += cb[l][inds[l] * vdim + j];
		}
	}
}
/****************************************************************************/
static void msvq_inv( float *y, int *inds, const MSVQ *msvq_tab )
/*  <- y			: reconstruction vector */
/*  -> inds		: array of indices */
/*  -> msvq_tab	: all required parameters and tables */
{
	int           j, l;
	const int     vdim = msvq_tab->vdim;
	const int     stages = msvq_tab->nstages;
	const float **cb = msvq_tab->cbs;
	for( j = 0; j < vdim; j++ ) {
		y[j] = cb[0][inds[0] * vdim + j];
		for( l = 1; l < stages; l++ ) {
			y[j] += cb[l][inds[l] * vdim + j];
		}
	}
}

/*
	Predictive ms vector quantizer
*/
void pmsvq( float *y, int **prm, float *x, float *old_x, const PMSVQ *filt_hi_pmsvq )
{
	float        e[MAX_VECT_DIM];
	float        eq[MAX_VECT_DIM];
	int          i;
	const float *cbm = filt_hi_pmsvq->mean;
	float        a = filt_hi_pmsvq->a;
	int          n = filt_hi_pmsvq->msvq.vdim;
	int *        inds = *prm;
	*prm += filt_hi_pmsvq->msvq.nstages;
	/* compute the predictor error */
	for( i = 0; i < n; i++ ) {
		e[i] = ( x[i] - cbm[i] ) - a * old_x[i];
	}
	/* quantize the prediction error */
	msvq( eq, inds, e, &filt_hi_pmsvq->msvq );
	/* save for next frame */
	for( i = 0; i < n; i++ ) {
		old_x[i] = eq[i] + a * old_x[i];
		y[i] = old_x[i] + cbm[i];
	}
}
void pmsvq_inv( float *y, int **prm, float *old_y, const int bfi, const PMSVQ *filt_hi_pmsvq )
{
	float        eq[MAX_VECT_DIM];
	int          i;
	int *        inds = *prm;
	const float *cbm = filt_hi_pmsvq->mean;
	int          n = filt_hi_pmsvq->msvq.vdim;
	float        a = ( bfi ? filt_hi_pmsvq->a_fe : filt_hi_pmsvq->a );
	*prm += filt_hi_pmsvq->msvq.nstages;
	if( !bfi ) {
		/* dequantize the prediction error */
		msvq_inv( eq, inds, &filt_hi_pmsvq->msvq );
	}
	else {
		set_zero( eq, n );
	}
	/* reconstruct and save for next frame */
	for( i = 0; i < n; i++ ) {
		old_y[i] = eq[i] + a * old_y[i];
		y[i] = old_y[i] + cbm[i];
	}
}

/*
	Lpc synthesis
*/
void syn_filt(
    float a[],     /* input : LP filter coefficients                     */
    int   m,       /* input : order of LP filter                         */
    float x[],     /* input : input signal                               */
    float y[],     /* output: output signal                              */
    int   l,       /* input : size of filtering                          */
    float mem[],   /* in/out: initial filter states                      */
    int   update_m /* input : update memory flag: 0 --> no memory update */
    )              /*                             1 --> update of memory */
{
	int   i, j;
	float s, *yy;
	float buf[BUF_SIZ];
	yy = &buf[0];
	/* copy initial filter states into synthesis buffer */
	for( i = 0; i < m; i++ ) {
		*yy++ = mem[i];
	}
	for( i = 0; i < l; i++ ) {
		s = x[i];
		for( j = 1; j <= m; j++ ) {
			s -= a[j] * yy[i - j];
		}
		yy[i] = s;
		y[i] = s;
	}
	/* Update memory if required */
	if( update_m ) {
		for( i = 0; i < m; i++ ) {
			mem[i] = yy[l - m + i];
		}
	}
	return;
}
/*
	Lpc residual
*/
void residu(
    float *a, /* input : LP filter coefficients           */
    int    m, /* input : order of LP filter               */
    float *x, /* input : input signal (usually speech)    */
    float *y, /* output: output signal (usually residual) */
    int    l  /* input : size of filtering                */
)
{
	float s;
	int   i, j;
	for( i = 0; i < l; i++ ) {
		s = x[i];
		for( j = 1; j <= m; j++ ) {
			s += a[j] * x[i - j];
		}
		y[i] = s;
	}
	return;
}

/*
	Constrained cholseky linear equation solver
*/
int cholsolc( float r[HI_FILT_ORDER][HI_FILT_ORDER],
    float           c[HI_FILT_ORDER],
    float           h[HI_FILT_ORDER],
    int             n )
{
	float p[HI_FILT_ORDER];
	int   i, j, k;
	float sum;
	float lambda;
	float v[HI_FILT_ORDER];
	/* cholesky decomposition */
	for( i = 0; i < n; i++ ) {
		for( j = i; j < n; j++ ) {
			for( sum = r[i][j], k = i - 1; k >= 0; k-- ) {
				sum -= r[i][k] * r[j][k];
			}
			if( i == j ) {
				if( sum < EPSILON ) {
					return ( 1 );
				}
				p[i] = (float)sqrt( sum );
			}
			else {
				r[j][i] = sum / p[i];
			}
		}
	}
	/* linear system solving */
	for( i = 0; i < n; i++ ) {
		for( sum = c[i], lambda = 1.0f, k = i - 1; k >= 0; k-- ) {
			sum -= r[i][k] * h[k];
			lambda -= r[i][k] * v[k];
		}
		h[i] = sum / p[i];
		v[i] = lambda / p[i];
	}
	for( i = n - 1; i >= 0; i-- ) {
		for( sum = h[i], lambda = v[i], k = i + 1; k < n; k++ ) {
			sum -= r[k][i] * h[k];
			lambda -= r[k][i] * v[k];
		}
		h[i] = sum / p[i];
		v[i] = lambda / p[i];
	}
	sum = 0.0f;
	lambda = 0.0f;
	for( i = 0; i < n; i++ ) {
		sum += h[i];
		lambda += v[i];
	}
	lambda = -sum / lambda;

	sum = 0;
	for( i = 0; i < n; i++ ) {
		sum += h[i];
	}
	sum = (float)fabs( sum );
	for( i = 0; i < n; i++ ) {
		h[i] += lambda * v[i];
	}
	if( sum > 1.0 ) {
		sum = 1.0f / sum;
		for( i = 0; i < n; i++ ) {
			h[i] *= sum;
		}
	}
	return ( 0 );
}

float glev_s( float *b, /* output: filter coefficients */
    float *          r, /* input : vector of autocorrelations   */
    float *          z, /* input: vector of cross correlations */
    int              m  /* input : order of filter            */
)
{
	float  buf[( WIENER_ORDER + 1 ) * ( WIENER_ORDER + 1 )];
	float  a[( WIENER_ORDER + 1 ) * ( WIENER_ORDER + 1 )];
	float *rc; /* reflection coefficients  0,...,m-1 */
	float  s, at, err, t;
	int    i, j, l;
	rc = &buf[0];
	if( r[0] == 0 ) {
		r[0] = 1.0;
	}
	rc[0] = ( -r[1] ) / r[0];
	a[0] = 1.0;
	a[1] = rc[0];
	err = r[0] + r[1] * rc[0];
	b[0] = z[0] / r[0];
	for( i = 2; i <= m; i++ ) {
		s = 0.0;
		for( j = 0; j < i; j++ ) {
			s += r[i - j] * a[j];
		}
		t = 0.0f;
		for( j = 0; j < i - 1; j++ ) {
			t += r[i - 1 - j] * b[j];
		}
		rc[i - 1] = ( -s ) / ( err );
		b[i - 1] = ( z[i - 1] - t ) / err;
		for( j = 0; j < i - 1; j++ ) {
			b[j] += b[i - 1] * a[i - 1 - j];
		}
		for( j = 1; j <= ( i / 2 ); j++ ) {
			l = i - j;
			at = a[j] + rc[i - 1] * a[l];
			a[l] += rc[i - 1] * a[j];
			a[j] = at;
		}
		a[i] = rc[i - 1];
		err += rc[i - 1] * s;
		if( err <= 0.0 ) {
			err = 0.01f;
		}
	}
	return err;
}

void crosscorr( float *vec1,   /* (i)	:	Input vector 1			*/
    float *            vec2,   /* (i)	:	Input vector 2			*/
    float *            result, /* (o)	:	Output result vector	*/
    int                length, /* (i)	:	Length of input vectors */
    int                minlag, /* (i)	:	Minimum lag				*/
    int                maxlag  /* (i)	:	Maximum lag				*/
)
{
	int   i, j, k;
	float t;
	k = 0;
	for( i = minlag; i < 0; i++ ) {
		t = 0;
		for( j = 0; j < length + i; j++ ) {
			t += vec1[j] * vec2[j - i];
		}
		result[k++] = t;
	}
	for( i = 0; i < maxlag; i++ ) {
		t = 0;
		for( j = 0; j < length - i; j++ ) {
			t += vec1[j + i] * vec2[j];
		}
		result[k++] = t;
	}
	return;
}
