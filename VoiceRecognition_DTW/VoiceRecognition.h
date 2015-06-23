#ifndef _VOICE_RECOGNITION_H_
#define _VOICE_RECOGNITION_H_
int detectEndPoint(float** f,int frameLen,int frameN,int *pBeginIndx,int *pEndIndex);
float DTW_Origin(float **A,int lenA,float **B,int lenB,float (*calDis)(float*,float*));
float calDis(float *a,float *b);
#define DTW(A,lenA,B,lenB) DTW_Origin(A,lenA,B,lenB,calDis)
#define PUTOUT_PATH 
float **WAVtoMFCCs(char *WavfileName,int *pLen);
void MFCCstoMFCCFile(float **MFCCs,int len,char * MFCCFileName);
void WAVtoMFCCFile(char *WavfileName);
float** MFCCfiletoMFCCs(char *MFCCFileName,int * pLen); 
void freeMFCCsMem(float **mfcc,int len);
extern char *gTemplateDir;
#endif