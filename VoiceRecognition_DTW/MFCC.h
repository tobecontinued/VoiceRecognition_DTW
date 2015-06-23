#ifndef _MFCC_H_
#define _MFCC_H_
#define TRIANGLE_Num (20)
#define DCT_Num (12)
#define MFCC_DIMENSION (39) //(12+1)*3
void preemphasize(float *sample, int sampleN);
float** framing(float *sample, int sampleN,int frameLen,int frameOffest,int* pframeN);
void freeFrames(float **f,int frameN);
void HammingWindowing(float *frame,int windowLen);
void FFT(float x[],float y[],int n,int sign);
float* MelScaleTriangleFilters(float *F,int n,int SampleRate,int k);
void DCT2(float *d,int n);
float* FrametoMFCC(float* f,int frameLen,int SampleRate);
#endif
