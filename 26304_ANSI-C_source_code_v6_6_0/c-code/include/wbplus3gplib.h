#ifndef _WBPLUS3GPLIB
#define _WBPLUS3GPLIB

int Create3GPAMRWBPlus( void );
int Create3GPAMRWB( void );
int WriteSamplesAMRWBPlus( EncoderConfig conf, void *Serial, int length );
int Close3GP( char *filename );
int GetNextFrame3GP( short *tfi, int *bfi, short *extension, short *mode, short *st_mode, short *fst, void *serial, int init );
// int GetNextFrame3GP(short *extension, short *mode, short *st_mode, short *fst,void *serial, int init); // old
int Open3GP( short *tfi, int *bfi, char *filename, int verbose, DecoderConfig *conf );
// int Open3GP(char *filename, int verbose, DecoderConfig *conf); // old
#endif
