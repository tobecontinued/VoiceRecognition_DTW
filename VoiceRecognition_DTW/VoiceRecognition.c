#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <float.h>//FLT_MAX��
#include  <string.h>
#include "VoiceRecognition.h"
#include "MFCC.h"
#include "AudioData.h"
static int calZeroCrossingRate(float *f,int frameLen);
static float calFrameEnergy(float *f,int frameLen);
/*
 * ˫���޷�
 */
static  float MAX(float a,float b)
{
	return a > b ? a : b;
}
static float MIN(float a,float b)
{
	return a < b ? a : b;
}
static int ROUND(float x)
{
	return (int)((x) < 0 ? (x)-0.5 : (x)+0.5);
}
int detectEndPoint(float** f,int frameLen,int frameN,int *pStartIndx,int *pEndIndex)
{
	float *FrameEnergy,E_Tl,E_Th,Emax,Emin,Eave,Zave;
	int i,*ZCR,ZCR_Tl,ZCR_Th,Zmax,
		status,// %״̬��0������,���ɶ�,1������,2������
		maxSilence = 8,//�������϶ʱ��
		minAudio = 16,//�������ʱ��      
		holdTime = 0,//��������ʱ��
		silenceTime = 0;//������϶ʱ��
	ZCR = (int*)calloc(frameN,sizeof(int));
	FrameEnergy = (float*)calloc(frameN,sizeof(float));
	Zave = Zmax = ZCR[0] =  calZeroCrossingRate(f[0],frameLen);
	Emax = Emin = Eave = FrameEnergy[0] =  calFrameEnergy(f[0],frameLen);
	for(i = 1;i < frameN;++i)
	{
		ZCR[i] = calZeroCrossingRate(f[i],frameLen);
		Zave += ZCR[i];
		if(Zmax  <  ZCR[i] )
			Zmax = ZCR[i];
		FrameEnergy[i] = calFrameEnergy(f[i],frameLen);
		Eave += FrameEnergy[i];
		if(Emax  <  FrameEnergy[i])
			Emax = FrameEnergy[i];
		if(Emin > FrameEnergy[i] )
			Emin = FrameEnergy[i];
	}
	Zave /= frameN;
	Eave /= frameN;
	ZCR_Th	= MAX( ROUND(Zmax * 0.1) , 9);//5
	ZCR_Tl = MAX( ROUND(Zave * 0.1 ) , 3);//3
	E_Th = (float)MAX(Emin * 10, MAX(Eave * 0.2, Emax * 0.1));
	E_Tl = (float)MIN(Emin * 10, MIN(Eave * 0.2, Emax * 0.1));
	status = 0;
	for(i = 0;i < frameN && status !=2;++i)
	{
		switch(status)
		{
			//������.���ȶ�		
			case 0:
				if(ZCR[i] > ZCR_Th || FrameEnergy[i] > E_Th)
				{
					status = 1;
					*pStartIndx = i - holdTime;
					++holdTime;
					silenceTime = 0;
				}
				else if(ZCR[i] > ZCR_Tl || FrameEnergy[i] > E_Tl)
				{
					++holdTime;
				}
				else
				{
					holdTime = 0;
				}
				break;
			//������	
			case 1:
				if(ZCR[i] > ZCR_Tl || FrameEnergy[i] > E_Tl)
				{
					++holdTime;
				}
				else
				{
					++silenceTime;
					if(silenceTime < maxSilence)
					{
						++holdTime;
					}
					else if((holdTime-silenceTime) < minAudio)
					{
						status = 0;
						holdTime = 0;
						silenceTime = 0;
					}
					else
					{
						status = 2;
					}
				}
				break;
		}
	}
	holdTime -= silenceTime;
	*pEndIndex = *pStartIndx + holdTime;
	if (*pEndIndex == frameN) --*pEndIndex;
	free(ZCR);
	free(FrameEnergy);
	if(!holdTime)
		return 0;
	return 1;
}
float calDis(float *a,float *b)
{
	float dis = 0;
	int i;
	for(i = 0;i <  MFCC_DIMENSION;++i)
	{
		float diff = a[i];
		diff -=b[i];
		dis += diff * diff;
	}
	return (float)sqrt(dis);
}
/*
 * ��D(i,j)Ϊ�ۻ����룬d(i.j)Ϊ�㣨i��j���ľ���
 * ��̬����ĵݹ�ʽ D(n,m) = d(n,m) + min(D(n-1,m),D(n-1,m-1),D(n-1,m-2))
							 d(n,m) (n=0,m=0,1,2 �� m=0,n=0,1,2)
 * input : A   :ģ����������A[lenA][Dimension]
 *         lenA:A�ĳ���
		   B   :������������ A[lenA][Dimension]
		   lenB��
		   calDis:ָ�����2�������ľ���ĺ���
	output:�ۻ����� D(lenA-1,lenB-1)
 */
float DTW_Origin(float **A,int lenA,float **B,int lenB,float (*calDis)(float*,float*))
{
	//�þ��� D = (lenA + 1 ) * (lenB + 1 )
	//�������ϵ���
	 //�ݹ�ת���ɵ���
	int i,j;
	float **D,d;
	D = (float**)calloc(lenA + 1,sizeof(float*));
	for(i = 0;i <= lenA;++i)
	{
		D[i] = (float*)calloc(lenB + 1,sizeof(float)); 
		for(j = 0;j <= lenB;++j)
		{
			if (i == 0 || j == 0)
				D[i][j] = FLT_MAX;
			else
				D[i][j] = calDis(A[i - 1],B[j - 1]);
		}
	}
	D[0][0] = 0.0;

	for(i = 0;i  <  lenA; ++i)
	{
		for (j = 0; j < lenB; ++j)
		{
				d = MIN( MIN(D[i][j], D[i][j + 1]), D[i + 1][j]);
				if (d == FLT_MAX)
					D[i + 1][j + 1] = FLT_MAX;
				else
					D[i + 1][j + 1] += d;
		}
	}
	

	d = D[lenA][lenB];
	for(i = 0;i <= lenA;++i)
	{
		free(D[i]);
	}
	free(D);
	
	return	d;
}
/* ����һ֡�Ĺ�����
 * ������ʽ
 * ZCR = 0.5*sum(sgn(f[i]) -sgn(f[i-1]),i=1,n-1)
 * ʵ��һ�ι��㶨��Ϊ��xn(m)*xn(m-1) < 0 && |xn(m) - xn(m-1)| > ��
 * 
 */
#define SGN(x) (x>=0?1:-1)
int calZeroCrossingRate(float *f,int frameLen)
{
	int ZRC = 0,i;
	float theta = 0.01f;
	for(i = 1;i < frameLen;++i)
	{
		if((SGN(f[i]) - SGN(f[i-1])) && fabs(f[i] - f[i -1]) > theta)
			++ZRC;
	}
	return ZRC / 2;
}
/*
 * En = sum(f[i]^2,i=0,n-1)
 */
float calFrameEnergy(float *f,int frameLen)
{
	float E = 0;
	int i;
	for(i = 0;i < frameLen;++i)
	{
		E += f[i] * f[i];
	}
	return E;
}
/*
 * ��ȡwav�ļ������MFCC
 * input :WavfileName :��Ƶ�ļ���
 * output:len MFCC�ĳ���
 * return :MFCC���� [len][dismension]
 */
float **WAVtoMFCCs(char *WavfileName,int *pLen)
{
	WavFileHeaderStruct *fileHeader ;
	float* WAVData,**f,**mfcc;
	int i,n,frameLen,frameN,start,end;
	int sampleRate;
	fileHeader = openWAVFile(WavfileName);
	if (!fileHeader)
	{
		fprintf(stderr, "error��%s:%d\n", __FILE__, __LINE__);
		return NULL;
	}
	WAVData = readAllWAVData(fileHeader, &n);
	sampleRate = fileHeader->FormatChunk.SampleRate;
	closeWAVFIle(fileHeader);
	
	//Ԥ��ǿ
	preemphasize(WAVData,n);
    //��֡
	if(sampleRate  < 16000)
	{
		frameLen = 256;
	}
	else 
	{//16000
		frameLen = 512;
	}
	f = framing(WAVData,n,frameLen,frameLen/3,&frameN);
	free(WAVData);
	
	//��˵㣺��ÿ֡��mfccϵ��
    if(detectEndPoint(f,frameLen,frameN,&start,&end))
	{
		mfcc = (float**)calloc(end - start + 1,sizeof(float*));
		for(i = start;i  <= end;++i)
		{
			mfcc[i - start] = FrametoMFCC(f[i],frameLen,sampleRate);
		}
		*pLen = end - start + 1;
		return mfcc;
	}
	fprintf(stderr, "error��%s:%d\n", __FILE__, __LINE__);
	return NULL;
}
#ifdef WIN32
char *gTemplateDir = "template\\";
#else
char *gTemplateDir = "template/";
#endif
/*
 * ����ȡ��ģ����Ƶ��MFCC��������������ļ�
 * ������ gTemplateDir��
 * ģ���ļ���ʽ�� (int)len + float(mcff[0][0]...mcff[0][dimension-1],mcff[1][0]..mcff[len-1][dimension-1]) 
 */
void MFCCstoMFCCFile(float **MFCCs,int len,char * MFCCFileName)
{
	FILE *pfile;
	char path[64] ;
	strcpy(path,gTemplateDir);
	strcat(path,MFCCFileName);
	if(pfile = fopen(path,"wb"))
	{
		int i;
		fwrite(&len,sizeof(float),1,pfile);
		for(i = 0;i < len;++i)
		{
			fwrite(MFCCs[i],sizeof(float),MFCC_DIMENSION,pfile);
		}
		fclose(pfile);
	}
	else
	{
		//���ļ�ʧ��
		fprintf(stderr, "error��%s:%d\n", __FILE__, __LINE__);
	}
}
void WAVtoMFCCFile(char *WavfileName)
{
	int len,i;
	float ** MFCCs;
	char *MFCCsFileName = (char*)calloc(strlen(WavfileName) + 1,sizeof(char));
	i = strlen(WavfileName) - 1;
#ifdef WIN32
	while( i >=0 && WavfileName[i--] != '\\');
#else
	while( i >=0 && WavfileName[i--] != '/');
#endif
	i = i < 0?0:i;
	strcpy(MFCCsFileName,WavfileName + i + 2);
	i = strlen(MFCCsFileName) - 1;
	while( i >=0 && MFCCsFileName[i--] != '.');
	if(i > 0)
	{
		++i;
		MFCCsFileName[++i] = 'm';
		MFCCsFileName[++i] = 'f';
		MFCCsFileName[++i] = 'c';
		MFCCs = WAVtoMFCCs(WavfileName,&len);
		MFCCstoMFCCFile(MFCCs,len,MFCCsFileName);
	}
	else
	{
		//�ļ������淶
		fprintf(stderr, "error��%s:%d\n", __FILE__, __LINE__);
	}
	free(MFCCsFileName);
}
/*
 * ��ȡģ���ļ�,������Ӧ��MFCCs����
 */
float** MFCCfiletoMFCCs(char *MFCCFileName,int * pLen)
{
	float **MFCCs = NULL;
	FILE *pfile;
	if(pfile =fopen(MFCCFileName,"rb"))
	{
		int i;
		fread(pLen,sizeof(int),1,pfile);
		MFCCs = (float**)calloc(*pLen,sizeof(float*));
		for(i = 0;i < *pLen;++i)
		{
			MFCCs[i] = (float*)calloc(MFCC_DIMENSION,sizeof(float));
			fread(MFCCs[i],sizeof(float),MFCC_DIMENSION,pfile);
		}
	}
	else
	{
	}
	return MFCCs;
}
/*
 * �ͷ�mfcc�Ķ�̬�ڴ�
 */
void freeMFCCsMem(float **mfcc,int len)
{
	int i;
	for(i = 0;i < len;++i)
		free(mfcc[i]);
	free(mfcc);
}
