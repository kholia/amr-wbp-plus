#ifndef _NCLASS_H_
#define _NCLASS_H_

#include "../lib_amr/typedef.h"

#define TWINLEN ( 24 - 8 )
#define LPHLEN 4
#define LPHAVELEN 8

#define TWINLENSHORT ( 4 )
#define COMPLEN 12

// For debugging
#define COMPLEN2 12
#define CS 0

#define ACELP_MODE 0
#define TCX_MODE 1
#define TCX_OR_ACELP 2

typedef struct
{

	float  levelHist[TWINLEN][COMPLEN];
	float  averageHistTime[COMPLEN];
	float  stdDevHistTime[COMPLEN];
	float  averageHistTimeShort[COMPLEN];
	float  stdDevHistTimeShort[COMPLEN];
	float  lphBuf[LPHLEN];
	float  lphAveBuf[LPHAVELEN];
	short  prevModes[4];
	Word16 vadFlag[4];
	Word16 vadFlag_old[4];
	Word16 LTPLag[10];
	float  NormCorr[10];
	float  LTPGain[10];
	float  TotalEnergy[5];
	short  NoMtcx[2];
	short  NbOfAcelps;
	float  ApBuf[4 * M];
	float  lph[4];
	short  StatClassCount;

	Word16 LTPlagV[8];

} NCLASSDATA;

short classifyExcitation( NCLASSDATA *stClass, float level[], short sfIndex );
void  classifyExcitationRef( NCLASSDATA *stClass, float *ISPs, short *coding_mod );

void initClassifyExcitation( NCLASSDATA *stClass );

#endif /*_NCLASS_H_*/
