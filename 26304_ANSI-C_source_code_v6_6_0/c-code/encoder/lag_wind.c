/*-------------------------------------------------------------*
 * procedure init_lag_wind:                                    *
 *           ~~~~~~~~~~~~~~                                    *
 * Compute the lag window (for windowing the autocorrelations  *
 * in LP analysis).                                            *
 *-------------------------------------------------------------*/
#include "../include/amr_plus.h"
/*-------------------------------------------------------------*
 * procedure lag_wind:                                         *
 *           ~~~~~~~~~                                         *
 * lag windowing of the autocorrelations                       *
 *-------------------------------------------------------------*/
void lag_wind(
    float r[], /* in/out: autocorrelations            */
    int   m    /* input : order of LP filter          */
)
{
	int i;
	for( i = 0; i <= m; i++ ) {
		r[i] *= lag_window[i];
	}
	return;
}
