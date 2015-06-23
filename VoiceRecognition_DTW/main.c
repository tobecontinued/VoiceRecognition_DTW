/*
 * 基于DTW的孤词识别
 * author: tbc.dengwenqi@gmail.com
 * all right reserved
 */
//听觉感知系统的特效
//http://i-math.sysu.edu.cn/Multimedia/multi/course1-9-1.html
//Audio Signal Processing and Recognition (音訊處理與辨識)
//http://neural.cs.nthu.edu.tw/jang/books/audioSignalProcessing/	
//Data Clustering and Pattern Recognition (資料分群與樣式辨認)
//http://neural.cs.nthu.edu.tw/jang/books/dcpr/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "VoiceRecognition.h"
#include "AudioData.h"
#include "MFCC.h"
//50.wav的MFCC求取，程序报错
//成功率不高
//三角滤波器组Ok，MFCC的改进 2015.4.4 22:38
//修改完成，80%正确： 2015.5.19 12;30
void storeTemplate()
{
	int i;
	char nameBuff[32];
	for (i = 0; i < 10; ++i)
	{
#ifdef WIN32
		sprintf(nameBuff,"train\\%d0.wav",i);
#else
		sprintf(nameBuff,"train/%d0.wav",i);
#endif
		WAVtoMFCCFile(nameBuff);
	}
}
int main()
{
	int i,j,k;
	char nameBuff[32];
	float **T, **t, min, dtwVal;
	int Tlen, tlen, minIndex, right = 0;
	storeTemplate();
	for (i = 0; i < 10; ++i)
	{
#ifdef WIN32
		sprintf(nameBuff,"test\\%d1.wav",i);
#else
		sprintf(nameBuff,"test/%d1.wav",i);
#endif

		t = WAVtoMFCCs(nameBuff, &tlen);
		minIndex = -1;
		min = FLT_MAX;
		if (t == NULL)
		{
			fprintf(stderr, "error：%s:%d\n", __FILE__, __LINE__);
			exit(-1);
		}
		for (j = 0; j < 10; ++j)
		{
#ifdef WIN32
		sprintf(nameBuff, "template\\%d0.mfc", j);
#else
		sprintf(nameBuff, "template/%d0.mfc", j);
#endif
			T = MFCCfiletoMFCCs(nameBuff, &Tlen);
			if (T == NULL)
			{
			//	fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
				exit(-1);
			}

			dtwVal = DTW(T, Tlen, t, tlen);
			if (min > dtwVal)
			{
				min = dtwVal;
				minIndex = j;
			}
			//		printf("test:%d train:%d dist:%f\n", i, j, dtwVal);
			freeMFCCsMem(T, Tlen);
		}
		printf("%d match %d,min dist:%f\n", i, minIndex, min);
		if (i == minIndex) ++right;
		freeMFCCsMem(t, tlen);
	}
	printf("right:%d , %%%f\n",right,(float)right / 10   * 100.0 );
	return 0;
}
