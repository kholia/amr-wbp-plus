/*---------------------------------------------------------------------*
 * routine int_lpc()                                                   *
 * ~~~~~~~~~~~~~~~~~                                                   *
 * Find the interpolated LPC parameters in every subframes.            *
 *---------------------------------------------------------------------*/
#include "../include/amr_plus.h"
void int_lpc_np1(
    float isf_old[], /* input : LSFs from past frame              */
    float isf_new[], /* input : LSFs from present frame           */
    float a[],       /* output: LP coefficients in both subframes */
    int   nb_subfr,  /* input: number of subframe                 */
    int   m          /* input : order of LP filter                */
)
{
	float isf[M], inc, fnew, fold, *p_a;
	int   i, k;
	inc = 1.0f / (float)nb_subfr;
	p_a = a;
	fnew = 0.0f;
	for( k = 0; k < nb_subfr; k++ ) {
		fold = 1.0f - fnew;
		for( i = 0; i < m; i++ ) {
			isf[i] = isf_old[i] * fold + isf_new[i] * fnew;
		}
		fnew += inc;
		E_LPC_f_isp_a_conversion( isf, p_a, m );
		p_a += ( m + 1 );
	}
	/* estimated coef for next subframe: use isf_new */
	E_LPC_f_isp_a_conversion( isf_new, p_a, m );
	return;
}
