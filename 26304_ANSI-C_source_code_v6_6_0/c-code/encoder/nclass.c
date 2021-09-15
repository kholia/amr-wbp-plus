#include "../include/amr_plus.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LOWLIMIT 8
#define LOWOMIT 1
void initClassifyExcitation( NCLASSDATA *stClass )
{
	int i, j;
	stClass->prevModes[0] = ACELP_MODE;
	stClass->prevModes[1] = ACELP_MODE;
	stClass->prevModes[2] = ACELP_MODE;
	stClass->prevModes[3] = ACELP_MODE;
	stClass->vadFlag[0] = 1;
	stClass->vadFlag[1] = 1;
	stClass->vadFlag[2] = 1;
	stClass->vadFlag[3] = 1;
	stClass->vadFlag_old[0] = 0;
	stClass->vadFlag_old[1] = 0;
	stClass->vadFlag_old[2] = 0;
	stClass->vadFlag_old[3] = 0;
	for( i = 0; i < 10; i++ ) {
		stClass->LTPGain[i] = 0;
		stClass->LTPLag[i] = 0;
		stClass->NormCorr[i] = 0;
	}
	stClass->TotalEnergy[0] = 0;
	stClass->TotalEnergy[1] = 0;
	stClass->TotalEnergy[2] = 0;
	stClass->TotalEnergy[3] = 0;
	stClass->TotalEnergy[4] = 0;
	stClass->NoMtcx[0] = 0;
	stClass->NoMtcx[1] = 0;
	stClass->NbOfAcelps = 0;
	stClass->lph[0] = 0;
	stClass->lph[1] = 0;
	stClass->lph[2] = 0;
	stClass->lph[3] = 0;
	for( j = 0; j < 4 * M; j++ ) {
		stClass->ApBuf[j] = 0.0;
	}
	for( j = 0; j < TWINLEN; j++ ) {
		for( i = 0; i < COMPLEN; i++ ) {
			stClass->levelHist[j][i] = 0.000001f;
		}
	}
	for( i = 0; i < COMPLEN; i++ ) {
		stClass->averageHistTime[i] = 1.0f;
		stClass->stdDevHistTime[i] = 1.0f;
		stClass->averageHistTimeShort[i] = 1.0f;
		stClass->stdDevHistTimeShort[i] = 1.0f;
	}
	for( i = 0; i < LPHLEN; i++ )
		stClass->lphBuf[i] = 6.0f;
	for( i = 0; i < LPHAVELEN; i++ )
		stClass->lphAveBuf[i] = 6.0f;
	stClass->StatClassCount = 0;
}
short classifyExcitation( NCLASSDATA *stClass, float levelNew[], Word16 sfIndex )
{
	float levH = 0.0, levL = 0.0;
	float tmp1, tmp2, stdAve, stdAveShort, lph;
	int   i, j;
	short headMode;
	float NormalizedBands[COMPLEN2] = { 0.0 }, bckr_est[COMPLEN2] = { 0.0 };
	float AverFreqFiltBand = 0.0;
	for( i = CS; i < COMPLEN2; i++ ) {
		levelNew[i] = (float)( levelNew[i] + 0.00001 );
		levelNew[i] = (float)( levelNew[i] / bw[i] );
		if( stClass->vadFlag[sfIndex] ) {
			stClass->levelHist[0][i] = levelNew[i];
		}
	}
	if( stClass->vadFlag[sfIndex] == 0 ) {
	}
	else {
		levH = 0.0f;
		levL = 0.0f;
		for( i = CS; i < COMPLEN2; i++ ) {
			stClass->averageHistTime[i] = stClass->averageHistTime[i] - stClass->levelHist[TWINLEN - 1][i] / TWINLEN;
			stClass->averageHistTime[i] = stClass->averageHistTime[i] + stClass->levelHist[0][i] / TWINLEN;
			stClass->averageHistTimeShort[i] = stClass->averageHistTimeShort[i] - stClass->levelHist[TWINLENSHORT - 1][i] / TWINLENSHORT;
			stClass->averageHistTimeShort[i] = stClass->averageHistTimeShort[i] + stClass->levelHist[0][i] / TWINLENSHORT;
			tmp1 = 0.0f;
			tmp2 = 0.0f;
			for( j = 0; j < TWINLEN; j++ ) {
				tmp1 = ( stClass->averageHistTime[i] - stClass->levelHist[j][i] );
				tmp2 = tmp2 + tmp1 * tmp1;
			}
			tmp2 = (float)sqrt( tmp2 / TWINLEN );
			tmp2 = (float)( tmp2 / stClass->averageHistTime[i] );
			stClass->stdDevHistTime[i] = tmp2;
			tmp1 = 0.0f;
			tmp2 = 0.0f;
			for( j = 0; j < TWINLENSHORT; j++ ) {
				tmp1 = ( stClass->averageHistTimeShort[i] - stClass->levelHist[j][i] );
				tmp2 = tmp2 + tmp1 * tmp1;
			}
			tmp2 = (float)sqrt( tmp2 / TWINLENSHORT );
			tmp2 = (float)( tmp2 / stClass->averageHistTimeShort[i] );
			stClass->stdDevHistTimeShort[i] = tmp2;
			if( i < LOWOMIT )
				;
			else if( i < LOWLIMIT ) {
				levL = levL + stClass->levelHist[0][i];
			}
			else {
				levH = levH + stClass->levelHist[0][i];
			}
		}
	}
	/*Analysis of energy & spectral variation*/
	for( i = 4; i > 0; i-- ) {
		stClass->TotalEnergy[i] = stClass->TotalEnergy[i - 1];
	}
	stClass->TotalEnergy[0] = 0.0;
	for( i = CS; i < COMPLEN2; i++ ) {
		if( levelNew[i] < bckr_est[i] ) {
			bckr_est[i] = levelNew[i];
		}
		stClass->TotalEnergy[0] += ( levelNew[i] - bckr_est[i] );
	}
	for( i = CS; i < COMPLEN2; i++ ) {
		NormalizedBands[i] = (float)( ( levelNew[i] - bckr_est[i] ) / ( stClass->TotalEnergy[0] + 0.0001 ) );
	}
	/*the mean frequency along filter bands */
	for( i = CS; i < COMPLEN2; i++ ) {
		AverFreqFiltBand += (float)( NormalizedBands[i] * VadFiltBandFreqs[i] );
	}
	levL = levL / 2200.0f;
	levH = levH / 4000.0f;
	stdAve = 0.0f;
	stdAveShort = 0.0f;
	for( j = CS; j < COMPLEN2; j++ ) {
		stdAve += stClass->stdDevHistTime[j];
		stdAveShort += stClass->stdDevHistTimeShort[j];
	}
	stdAve = stdAve / (float)( COMPLEN2 - CS );
	stdAveShort = stdAveShort / (float)( COMPLEN2 - CS );
	if( stClass->vadFlag[sfIndex] ) {
		for( i = CS; i < COMPLEN2; i++ ) {
			for( j = TWINLEN - 1; j > 0; j-- ) {
				stClass->levelHist[j][i] = stClass->levelHist[j - 1][i];
			}
		}
	}
	if( stClass->vadFlag[sfIndex] ) {
		for( j = LPHLEN - 1; j > 0; j-- ) {
			stClass->lphBuf[j] = stClass->lphBuf[j - 1];
		}
		stClass->lphBuf[0] = levL / levH;
	}
	lph = 0.0f;
	for( j = 0; j < LPHLEN; j++ ) {
		lph = lph + stClass->lphBuf[j];
	}
	lph = lph / (float)LPHLEN;
	if( ( stClass->vadFlag[sfIndex] ) && ( sfIndex == 3 ) ) {
		for( j = LPHAVELEN - 1; j > 0; j-- ) {
			stClass->lphAveBuf[j] = stClass->lphAveBuf[j - 1];
		}
	}
	stClass->lphAveBuf[0] = lph;
	lph = 0.0;
	for( j = 0; j < LPHAVELEN; j++ ) {
		lph = lph + stClass->lphAveBuf[j] * lwg[j];
	}
	/*Classifies the signal into three categories
	 based on the standard deviation and tilt of the signal*/
	/* Encode with ACELP = ACELP_MODE*/
	/* Encode with TCX   = TCX_MODE*/
	/* Encode with ACELP or TCX = TCX_OR_ACELP (decided in ClassB refinement)*/
	/*TCX*/
	if( stClass->StatClassCount == 0 ) {
		headMode = TCX_MODE;
		if( stdAve < 0.4 ) {
			headMode = TCX_MODE;
		}
		else if( lph > 280 ) {
			headMode = TCX_MODE;
		}
		else if( stdAve >= 0.4 ) {
			if( ( 5 + ( 1 / ( stdAve - 0.4 ) ) ) > lph ) {
				headMode = TCX_MODE;
			}
			else if( ( -90 * stdAve + 120 ) < lph ) {
				headMode = ACELP_MODE;
			}
			else {
				headMode = TCX_OR_ACELP;
			}
		}
	}
	else {
		headMode = TCX_OR_ACELP;
	}
	/*Noise like signal with TCX*/
	if( ( headMode == ACELP_MODE || headMode == TCX_OR_ACELP ) && AverFreqFiltBand > 2000 ) {
		headMode = TCX_MODE;
	}
	if( stClass->StatClassCount < 5 ) {
		if( headMode == TCX_OR_ACELP ) {
			if( stdAveShort < 0.2 ) {
				headMode = TCX_MODE;
			}
			else if( stdAveShort >= 0.2 ) {
				/* TCX */
				if( ( 2.5 + ( 1 / ( stdAveShort - 0.2 ) ) ) > lph ) {
					headMode = TCX_MODE;
				}
				/*ACELP*/
				else if( ( -90 * stdAveShort + 140 ) < lph ) {
					headMode = ACELP_MODE;
				}
				/*TCX_OR_ACELP*/
				else {
					headMode = TCX_OR_ACELP;
				}
			}
		}
	}
	if( headMode == TCX_OR_ACELP ) {
		if( stClass->StatClassCount < 15 ) {
			if( ( stClass->TotalEnergy[0] / stClass->TotalEnergy[1] ) > 25 ) {
				headMode = ACELP_MODE;
			}
		}
	}
	if( ( headMode == TCX_MODE || headMode == TCX_OR_ACELP ) ) {
		if( AverFreqFiltBand > 2000 && stClass->TotalEnergy[0] < 60 ) {
			headMode = ACELP_MODE;
		}
	}
	for( i = 2; i >= 0; i-- ) {
		stClass->lph[i + 1] = stClass->lph[i];
	}
	stClass->lph[0] = lph;
	if( stClass->vadFlag[sfIndex] == 0 ) {
		headMode = TCX_MODE;
	}
	return ( headMode );
}
void classifyExcitationRef( NCLASSDATA *stClass, float *ISPs, short *headMode )
{
	int   sfIndex, i, j;
	short SDminInd = 0, SDmaxInd = 0;
	int   tmp1 = 0, tmp2 = 0;
	float SD[4] = { 0.0 };
	float tmp1F;
#define DFTN 64
#define DFTNx2 128
#define LPC_N 16
	float ip[LPC_N] = { 0.0 }, mag[DFTN];
	float cos_t[DFTNx2], sin_t[DFTNx2];
	float x = 0, y = 0;
	for( sfIndex = 0; sfIndex < 4; sfIndex++ ) {
		for( i = 0; i < 4; i++ ) {
			SD[sfIndex] += (float)( fabs( ISPs[( sfIndex + 1 ) * M + i] - ISPs[sfIndex * M + i] ) );
		}
		if( SD[sfIndex] < SD[SDminInd] ) {
			SDminInd = (short)( sfIndex );
		}
		else if( SD[sfIndex] > SD[SDmaxInd] ) {
			SDmaxInd = (short)( sfIndex );
		}
	}
	/*In the case of switching, the history of buffers are updated with the values of current frame*/
	if( stClass->StatClassCount == 15 ) {
		sfIndex = 0;
		stClass->LTPLag[( sfIndex * 2 + 2 ) - 2] = stClass->LTPLag[( sfIndex * 2 + 2 )];
		stClass->LTPLag[( sfIndex * 2 + 2 ) - 1] = stClass->LTPLag[( sfIndex * 2 + 2 ) + 1];
		stClass->LTPGain[( sfIndex * 2 + 2 ) - 2] = stClass->LTPGain[( sfIndex * 2 + 2 )];
		stClass->LTPGain[( sfIndex * 2 + 2 ) - 1] = stClass->LTPGain[( sfIndex * 2 + 2 ) + 1];
		stClass->NormCorr[( sfIndex * 2 + 2 ) - 2] = stClass->NormCorr[( sfIndex * 2 + 2 )];
		stClass->NormCorr[( sfIndex * 2 + 2 ) - 1] = stClass->NormCorr[( sfIndex * 2 + 2 ) + 1];
	}
	for( sfIndex = 0; sfIndex < 4; sfIndex++ ) {
		if( stClass->vadFlag[sfIndex] != 0 && headMode[sfIndex] == TCX_OR_ACELP ) {
			if( SD[sfIndex] > 0.2 ) {
				headMode[sfIndex] = ACELP_MODE;
			}
			else {
				tmp1 = 0;
				tmp1 = abs( stClass->LTPLag[( sfIndex * 2 + 2 )] - stClass->LTPLag[( sfIndex * 2 + 2 ) + 1] );
				tmp2 = abs( stClass->LTPLag[( sfIndex * 2 + 2 ) - 2] - stClass->LTPLag[( sfIndex * 2 + 2 ) - 1] );
				if( tmp1 < 2 && tmp2 < 2 && abs( stClass->LTPLag[( sfIndex * 2 + 2 ) - 1] - stClass->LTPLag[( sfIndex * 2 + 2 )] ) < 2 ) {
					if( ( stClass->LTPLag[( sfIndex * 2 + 2 ) - 2] == 18 && stClass->LTPLag[( sfIndex * 2 + 2 ) + 1] == 18 ) || ( stClass->LTPLag[( sfIndex * 2 + 2 ) - 2] == 115 && stClass->LTPLag[( sfIndex * 2 + 2 ) + 1] == 115 ) ) {
						if( fabs( stClass->LTPGain[( sfIndex * 2 + 2 )] - stClass->NormCorr[( sfIndex * 2 + 2 )] ) < 0.1 && fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) + 1] - stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] ) < 0.1 && ( stClass->NormCorr[( sfIndex * 2 + 2 )] > 0.9 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] > 0.9 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) - 1] > 0.9 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) - 2] > 0.9 ) ) {
							headMode[sfIndex] = ACELP_MODE;
						}
						else {
							headMode[sfIndex] = TCX_MODE;
						}
					}
					else if( fabs( stClass->LTPGain[( sfIndex * 2 + 2 )] - stClass->NormCorr[( sfIndex * 2 + 2 )] ) < 0.1 && fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) + 1] - stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] ) < 0.1 && ( stClass->NormCorr[( sfIndex * 2 + 2 )] > 0.88 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] > 0.88 ) ) {
						headMode[sfIndex] = ACELP_MODE;
					}
					else if( fabs( stClass->LTPGain[( sfIndex * 2 + 2 )] - stClass->NormCorr[( sfIndex * 2 + 2 )] ) > 0.2 && fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) + 1] - stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] ) > 0.2 ) {
						headMode[sfIndex] = TCX_MODE;
					}
					else {
						if( sfIndex < 2 ) {
							stClass->NoMtcx[0]++;
						}
						else {
							stClass->NoMtcx[1]++;
						}
					}
				}
				if( headMode[sfIndex] == TCX_OR_ACELP ) {
					tmp1 = (int)( stClass->TotalEnergy[4] );
					for( i = 0; i < 4; i++ ) {
						if( (int)( stClass->TotalEnergy[i] ) > tmp1 ) {
							tmp1 = (int)( stClass->TotalEnergy[i] );
						}
					}
					if( tmp1 < 60 ) {
						if( SD[sfIndex] > 0.15 ) {
							headMode[sfIndex] = ACELP_MODE;
						}
						else {
							if( sfIndex < 2 ) {
								stClass->NoMtcx[0]++;
							}
							else {
								stClass->NoMtcx[1]++;
							}
						}
					}
				}
			}
		}
		else if( stClass->vadFlag[sfIndex] != 0 && headMode[sfIndex] == ACELP_MODE ) {
			if( ( abs( stClass->LTPLag[( sfIndex * 2 + 2 )] - stClass->LTPLag[( sfIndex * 2 + 2 ) + 1] ) < 2 ) && ( abs( stClass->LTPLag[( sfIndex * 2 + 2 ) - 2] - stClass->LTPLag[( sfIndex * 2 + 2 ) - 1] ) < 2 ) ) {
				if( ( stClass->NormCorr[( sfIndex * 2 + 2 )] < 0.80 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] < 0.80 ) && ( SD[SDmaxInd] < 0.1 ) ) {
					headMode[sfIndex] = TCX_MODE;
				}
			}
			if( ( stClass->lph[sfIndex] > 200 ) && ( SD[SDmaxInd] < 0.1 ) ) {
				headMode[sfIndex] = TCX_MODE;
			}
		}
	}
	for( sfIndex = 0; sfIndex < 4; sfIndex++ ) {
		/*80ms TCX is disabled if VAD is set to zero in last 80ms */
		if( ( stClass->vadFlag_old[sfIndex] == 0 ) && ( stClass->vadFlag[sfIndex] == 1 ) && ( headMode[sfIndex] == TCX_MODE ) ) {
			if( sfIndex < 2 ) {
				stClass->NoMtcx[0]++;
			}
			else {
				stClass->NoMtcx[1]++;
			}
		}
		if( headMode[sfIndex] != ACELP_MODE ) {
			tmp1F = 0;
			tmp1F += (float)fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) - 2] - stClass->NormCorr[( sfIndex * 2 + 2 ) - 2] );
			tmp1F += (float)fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) - 1] - stClass->NormCorr[( sfIndex * 2 + 2 ) - 1] );
			tmp1F += (float)fabs( stClass->LTPGain[( sfIndex * 2 + 2 )] - stClass->NormCorr[( sfIndex * 2 + 2 )] );
			tmp1F += (float)fabs( stClass->LTPGain[( sfIndex * 2 + 2 ) + 1] - stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] );
			tmp1F = (float)tmp1F / 4;
			if( tmp1F < 0.006 && ( stClass->NormCorr[( sfIndex * 2 + 2 )] > 0.92 ) && ( stClass->NormCorr[( sfIndex * 2 + 2 ) + 1] > 0.92 ) && ( stClass->LTPLag[( sfIndex * 2 + 2 )] > 21 ) && ( stClass->LTPLag[( sfIndex * 2 + 2 ) + 1] > 21 ) ) {
				for( i = 0; i < DFTNx2; i++ ) {
					cos_t[i] = (float)t_cos[i * N_MAX / DFTNx2];
					sin_t[i] = (float)t_sin[i * N_MAX / DFTNx2];
				}
				for( i = 0; i < LPC_N; i++ ) {
					ip[i] = stClass->ApBuf[sfIndex * M + i];
				}
				mag[0] = 0.0;
				for( i = 0; i < DFTN; i++ ) { /* calc DFT */
					x = y = 0;
					for( j = 0; j < LPC_N; j++ ) {
						x += ip[j] * cos_t[( i * j ) & ( DFTNx2 - 1 )];
						y += ip[j] * sin_t[( i * j ) & ( DFTNx2 - 1 )];
					}
					mag[i] = (float)( 1 / sqrt( x * x + y * y ) );
				}
				tmp1F = 0;
				for( i = 1; i < 40; i++ ) { /*First element left out*/
					tmp1F += mag[i];
				}
				if( tmp1F > 95 && mag[0] < 5 ) {
					headMode[sfIndex] = TCX_MODE;
				}
				else {
					headMode[sfIndex] = ACELP_MODE;
					if( sfIndex < 2 ) {
						stClass->NoMtcx[0]++;
					}
					else {
						stClass->NoMtcx[1]++;
					}
				}
			}
		}
		if( stClass->StatClassCount < 12 ) {

			if( headMode[sfIndex] == TCX_OR_ACELP ) {
				tmp1 = 0;
				tmp2 = 0;
				for( i = 0; i < 4; i++ ) {
					if( ( stClass->prevModes[i] == 3 || stClass->prevModes[i] == 2 ) && stClass->vadFlag_old[i] == 1 && stClass->TotalEnergy[i] > 60 ) {
						tmp1++;
					}
					if( stClass->prevModes[i] == ACELP_MODE ) {
						tmp2++;
					}
					if( sfIndex != i ) {
						if( headMode[i] == ACELP_MODE )
							tmp2++;
					}
				}
				if( tmp1 > 3 ) {
					headMode[sfIndex] = TCX_MODE;
				}
				else if( tmp2 > 1 ) {
					headMode[sfIndex] = ACELP_MODE;
				}
				else {
					headMode[sfIndex] = TCX_MODE;
				}
			}
		}
		else {

			if( headMode[sfIndex] == TCX_OR_ACELP ) {
				headMode[sfIndex] = TCX_MODE;
			}
		}
	}
	stClass->NbOfAcelps = 0;
	for( sfIndex = 0; sfIndex < 4; sfIndex++ ) {
		if( headMode[sfIndex] == ACELP_MODE ) {
			stClass->NbOfAcelps++;
		}
	}
	/*Buffer updates*/
	stClass->LTPGain[0] = stClass->LTPGain[8];
	stClass->LTPGain[1] = stClass->LTPGain[9];
	stClass->LTPLag[0] = stClass->LTPLag[8];
	stClass->LTPLag[1] = stClass->LTPLag[9];
	stClass->NormCorr[0] = stClass->NormCorr[8];
	stClass->NormCorr[1] = stClass->NormCorr[9];
	for( i = 0; i < 4; i++ ) {
		stClass->vadFlag_old[i] = stClass->vadFlag[i];
		if( stClass->StatClassCount > 0 ) {
			stClass->StatClassCount--;
		}
	}
	return;
}
