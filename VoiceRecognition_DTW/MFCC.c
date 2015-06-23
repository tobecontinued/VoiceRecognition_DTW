#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MFCC.h"
//MFCC主要参考来源
//http://baike.baidu.com/view/2930343.htm#2
//http://blog.csdn.net/c395565746c/article/details/6210920
//http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
//其中的三角滤波器特别重要
//http://blog.csdn.net/xiaoding133/article/details/8106672?reload

/*
 * 该过程可以达到在音框化阶段对静音数据的判断，
 * 因为静音数据的值是几乎不变的，所以在做差分以后值会很小，接近于0，而有声音的数据则会保留较大的值
 * 从频域的角度来看
 * 对于语音信号的频谱，通常是频率越高幅值越小，在语音信号的频率增加两倍时，其功率谱的幅度下降6dB
 * 因此必须对高频进行加重处理，一般是将语音信号通过一个一阶高通滤波器H(z) = 1-0.9375z^-1，即为预加重滤波器。其目的是滤除低频干扰
 * input : sample  ：信号的时间序列
           sampleN ：信号的时间序列的长度
 * output: sample  : 预计权重后的序列填充的数组，保证长度>=sampleN
 */
void preemphasize(float *sample, int sampleN)   
{     
   int i;   
   float emphFac = (float)0.97;    
   for (i = sampleN - 1; i > 0; --i) 
   {   
        sample[i] = sample[i] - emphFac * sample[i-1];
		//这个公式 Sr(k) = e(k)(1 - S^-1)（移序算子）与一阶有限激励响应高通滤波器 H(z) = 1-0.9375z^-1是等价的
		//H(s)  = 1- 0.9375e^-s （拉式变换）
		//H(jw) = 1- 0.9375e^-jw 
    }   
    sample[0] = (float)(1.0 - emphFac) * sample[0];
}
/*
 * 语音数据的分帧
 * 语音信号是一种典型的非平稳信号，它的均值函数u(x)和自相关函数R(xl,x2)都随时间而发生较大的变化。
 * 但研究发现，语音信号在短时间内频谱特性保持平稳，即具有短时平稳特性。
 * 因此，在实际处理时可以将语音信号分成很小的时间段(约10~30ms，一般取帧长20ms，帧移为帧长的1/3~1/2。)，称之为“帧”，
 * 作为语音信号处理的最小单位，帧与帧的非重叠部分称为帧移，而将语音信号分成若干帧的过程称为分帧
 * input : sample  :需要分帧的的序列
		   sampleN ：信号的时间序列的长度
		   frameLen ：每一帧的长度，最好为2的幂，方便FFT,DCT等计算
		   frameOffest : frameLen * 0.33
   output : pframeN ：帧数将填充此整数
   return : 该函数参数一个数组[帧数][帧长]，表示分帧后的数据,
 */
float** framing(float *sample, int sampleN,int frameLen,int frameOffest,int* pframeN)
{
	int i,j,k;
	float **framedSample;
	*pframeN = (sampleN - frameLen) / frameOffest;
    framedSample = (float**)calloc(*pframeN,sizeof(float*));
	for(i = 0, k =0; i < *pframeN;i++)
	{
		framedSample[i] = (float*)calloc(frameLen,sizeof(float));
		for(j = 0; j < frameLen;j++)
		{

			framedSample[i][j] = sample[k++];
		}
		k = k - frameLen + frameOffest;
	}
	return framedSample;
}
void freeFrames(float **f,int frameN)
{
	while(frameN)
	 free(f[--frameN]);
	free(f);
}
/*
 * 汉明窗函数，为了保持语音信号的短时平稳性，利用窗函数来减少由截断处理导致的Gibbs效应
 * 汉明窗函数:
      Whm(n) = 0.5-0.46*cos(2*PI*n/(N-1)) ((0≤n＜N-1))
	      = 0                             (Other)
 * Gibbs效应:对于有不连续点的函数,取傅里叶级数后，即使所取项数趋于无穷，在不连续处仍不收敛
			 与原函数，在越变点附近的波形，总是不可避免的存在有起伏振荡，从而在越变点附近某些
			 点的函数值超过1,而形成从过冲，随着级数所取项数的增多，起伏振荡存在的时间将缩短，
			 但过冲值趋于越9%的固定值--《信号与线性系统》，P104
 * 汉明窗函数的原理暂时不清楚
 * input  : frame :需要加窗的一帧数据
          ：windowLen 窗的长度，通常取帧长
 * output : frame 加窗后的数据
 *
 */
#define PI 3.1415
void HammingWindowing(float *frame,int windowLen)
{
	int i;
	for(i = 0;i < windowLen;i++)
	{
		frame[i] = frame[i] * (float)(0.5 - 0.46*cos(2*PI*i/(windowLen -1)));
	}
}
/*
 * 快速傅里叶变换 --《信号与线性系统》第九章
 * if(n == 2)
    F(0) = x(0) + x(1)
	F(1) = x(0) - x(1)
	return
 * for i =0; i < n/2;++i
 *  f1(i) = f(2i)
 *  f2(i) = f(2i+1)
 * F1(k) = FFFT(f1)
 * F2(k) = FFT(f2)
 * for i =0; i < n/2;++i
	F(i) = F1(i) + WnkF2(i)
	F(i+n/2) = F1(i) - WnkF2(i)
* 设fs为时序的抽样频率
*  FFT后Fo的频域的抽样间隔
* 则 Fo = fs/n	
	
 * input :  Real :时间序列的实部
            Imag :时间序列的虚部
		    n    :时间序列的长度,n应该满足 n == 2^k;
			不满足时补0值达到序列长度要求       
 * output : Real :傅里叶变换后的实部
            Imag :傅里叶变换后的虚部

* Wnk = e^(-j*2*pi*k/n)
*/
//数组x存储时域序列的实部，数组y存储时域序列的虚部
//n代表N点FFT，sign=1为FFT，sign=-1为IFFT
//http://my.csdn.net/lixeb/code/detail/32515
void FFT(float x[],float y[],int n,int sign)
{
	int i,j,k,l,m,n1,n2;
	float c,c1,e,s,s1,t,tr,ti;
	//Calculate i = log2N
	for(j = 1,i = 1; i < 16; i++)
	{
		 m = i;
		 j = 2*j;
		 if(j == n)
			 break;
	}
	//计算蝶形图的输入下标（码位倒读）
	n1 = n - 1;
	for(j=0,i=0; i < n1; i++)
	{
		if(i < j)           
		{
			 tr = x[j];
			 ti = y[j];
			 x[j] = x[i];
			 y[j] = y[i]; 
			 x[i] = tr;
			 y[i] = ti;                 
		}
		k = n/2;
		while(k < (j+1))
		{
			j = j - k;
			k = k/2;              
		}
		j = j + k;
	}
	//计算每一级的输出，l为某一级，i为同一级的不同群，使用同一内存（即位运算）
	n1 = 1;
	for(l=1; l <= m; l++)
	{
		n1 = 2*n1;
		n2 = n1/2;
		e = (float)(PI / n2);
		c = 1.0;
		s = 0.0;
		c1 = (float)cos(e);
		s1 = -sign*(float)sin(e);
		for(j=0; j < n2; j++)
		{
			for(i=j; i < n; i+=n1)         
			{
				k = i + n2;
				tr = c*x[k] - s*y[k];
				ti = c*y[k] + s*x[k];
				x[k] = x[i] - tr;
				y[k] = y[i] - ti;
				x[i] = x[i] + tr;
				y[i] = y[i] + ti;        
			}
			t = c;
			c = c*c1 - s*s1;
			s = t*s1 + s*c1;
		}
	}
	//如果是求IFFT，再除以N
	if(sign == -1)
	{
		for(i=0; i < n; i++)
		{
			x[i] /= n;
			y[i] /= n;
		}
	}
}
/*
 * 三角带通滤波器：在频域呈三角型
 * 设o,c,h分别为三角滤波器的下限，中心，上限频率y
 * on，cn，hn为频率除以频域采样间隔
 * 传函 W(k) = (k - on）/(cn-on) on < = k  < =cn
			   (hn-k) /(hn- cn)  cn < = k  < =hn	
 * input  : F	 :需要滤波的谱
            n	 :谱的长度
			Fbeta   : 频域相邻2点的频率差
		    oFreq :下限频率
		    hFreq :上限频率
 * output : v 
 */
void TriangleFilter(float *F,int Fbeta,int n,float oFreq,float hFreq)
{
	/*
	 * 设fs为时序的抽样频率
	 *  FFT后Fo的频于的抽样间隔
	 * 则 Fo = fs/n
	 */	
	int on,cn,hn,i;
	on = (int)oFreq / Fbeta;
	hn = (int)hFreq / Fbeta;
	cn = (on + hn) / 2;
	for(i = 0;i < on;++i)
		F[i] = 0;
	for(;i < cn;++i)
		F[i] *= (i - on) / (cn - on);
	for(;i < hn;++i)
		F[i] *= (hn - i) / (hn - cn);
	for(;i < n;++i)
		F[i] = 0;

}
/*
 * 三角滤波器覆盖的范围都近似于人耳的一个临界带宽,以此来 模拟人耳的掩蔽效应(具体来说应该是频域的掩蔽效应)
 * 对功率谱滤波
 * 设采样频率为fs，由抽样定理，信号的最高频率为fg为fs/2(用fg = fs/n *(n-1))
 * 有主观音高fmel = 2595 *log(1 + f /700hz)
 * f为频率
 * 得到最高的主观音高 fmax
 * 在Mel刻度范围内，各个三角滤波器的中心频率是相等间隔的线性分布
 * 可以计算两个相邻三角滤波器的中心频率的间距为
 * deta = fmax / (k+1)
 * 每一个三角形滤波器的中心频率c(l) __在Mel频率轴上等间隔分布__。设o(l),c(l),h(l) 分别是第l 个三角形滤波器的下限，中心，和上限频率，
 * 则相邻三角形滤波器之间的下限，中心，上限频率的关系如下：c(l)=h(l-1)=o(l+1)
 * input : F :需要滤波的功率谱
           n :长度
		   SampleRate：采样率
		   k :滤波器个数
 * return : 返回长度为K的的数组，每个值为通过每个滤波器的对数能量，定義為訊號的平方和，再取以 10 為底的對數值，再乘以 10
 * 输出应该为 k个参数
 * 第K参数为通过第K的滤波器后数据求和
 * 目前此函数不对
 */
float freqtoMel(float f)
{
	return 2596.0 * log10(1 + f / 700.0);
}
float meltoFreq(float m)
{
	return  700.0 * (pow(10, m / 2596.0) - 1);
}
float* MelScaleTriangleFilters(float *F,int n,int SampleRate,int k)
{
	float *d,
		  detaSamplefreq = (float)SampleRate / n ,
		  detaMel = SampleRate / 2  / (k + 1),
		  oMel    = 0,
		  cMel    = detaMel,
		  hMel    = detaMel * 2,
		  fMel = 0,f = 0;
	int i =0,j,midIndex;
	d = (float*)calloc(k,sizeof(float));
	for(j = 0;j < k && i  < n;++j)
	{
		d[j] = 0.0;
		while(fMel  <= hMel && i  < n)
		{
			if(fMel  <  cMel)
			{
				d[j] += F[i] *(fMel - oMel)/detaMel;
				midIndex = i;
			}
			else
			{
				d[j] += F[i] * (hMel - fMel) /detaMel;
			}
			++i;
			f += detaSamplefreq;
			fMel = (float)freqtoMel(f);
		}
		i = midIndex + 1 ;
		f = detaSamplefreq * (i-1);
		fMel = (float)freqtoMel(f);
		oMel += detaMel;
		cMel += detaMel;
		hMel += detaMel;
	}
	for(j = 0;j < TRIANGLE_Num;++j)
	{
		if(d[j])
			d[j] = (float)log10(d[j]) * 10;  
	}
	return d;

}
/*
 * 三角滤波器组
 * input: k : 滤波器组的个数
 *		  n : 需要滤波的序列长度
 *	sampleRate : 采样率
 * return ： 滤波器器组序列，为float[k][n]数组
*/
float** melbank(int k, int n, int sampleRate)
{
	int i,j;
	float melDelta = freqtoMel(sampleRate / 2.0) / (k + 1);//滤波器的中心频率间隔
	float freqDelta = sampleRate / n;
	float freq;
	float *f = (float*)calloc(k + 2,sizeof(float));//滤波器组的中心频率
	float **melbanks = (float**)calloc(k,sizeof(float*));
		
	for (i = 0; i < k; ++i)
	{
		melbanks[i] = calloc(n,sizeof(float));
		for (j = 0; j < n; ++j)
			melbanks[i][j] = 0.0;
	}
	//求滤波器组频域的中心频率，它在mel域为等距分布
	f[0] = 0.0;
	for (i = 1; i <= k + 1; ++i)
		f[i] = f[i - 1] + melDelta;
	for (i = 1; i <= k + 1; ++i)
		f[i] = meltoFreq(f[i]);
	//滤波器组
	for (i = 0; i < k; ++i)
	{
		freq = 0.0;
		for (j = 0; j < n / 2; ++j)
		{
			if (freq >= f[i] && freq < f[i + 1])
				melbanks[i][j] = (freq - f[i]) / (f[i + 1] - f[i]);
			else if (freq >= f[i + 1] && freq <= f[i + 2])
				melbanks[i][j] = (f[i + 2] - freq) / (f[i + 2] - f[i + 1]);
			freq += freqDelta;
		}
	}
	free(f);
	return melbanks;
}
/*
 * 求m[beg,end)的差分,放在其后
 */
void deltaCoeff(float *m, int beg, int end)
{
	int i,j;
	j = end;
	m[j++] = 0.0;
	m[j++] = 0.0;
	for (i = beg + 2; i < end - 2; ++i,++j)
	{
		m[j] = 2.0 * (m[i + 2] - m[i - 2]) + m[i + 1] - m[i - 1];
		m[j] /= sqrt(2.0 * (1.0 + 2.0* 2.0));
	}
	m[j++] = 0.0;
	m[j++] = 0.0;
}
/*
 * 每帧能量，再加12;为dct输出，组成13维MFCC
 * 在加上一次差分，和二次差分，共有39维
 */
float* FrametoMFCC(float* f,int frameLen,int SampleRate)
{
	int i,j;
	float *fImag,*mfcc,
		  logEnergy = 0.0,
		  *logFilterEnergy;
	float **melBanks;
	mfcc = (float*)calloc(MFCC_DIMENSION,sizeof(float));
	//加窗
	HammingWindowing(f ,frameLen);
	//计算一帧的对数能量 定義為一個音框內訊號的平方和，再取以 10 為底的對數值，再乘以 10
	for(j = 0;j < frameLen;++j)
	{
		logEnergy += f[j] * f[j];
	}

	logEnergy = (float)log10(logEnergy) * 10;
	//计算功率谱
	fImag = (float*)calloc(frameLen,sizeof(float));
	for(j = 0;j < frameLen;++j)
	{
		fImag[j] = 0.0;
	}
	FFT(f,fImag ,frameLen,1);
	for(j =0;j < frameLen;++j)
	{
		f[j] = f[j] * f[j] + fImag[j] * fImag[j];
	}
	free(fImag);
	//将功率谱通过滤波器组,并求每个滤波器的对数能量 
	logFilterEnergy = calloc(TRIANGLE_Num,sizeof(float));
	melBanks = melbank(TRIANGLE_Num,frameLen,SampleRate);
	for (i = 0; i < TRIANGLE_Num; ++i)
	{
		logFilterEnergy[i] = 0.0;
		for (j = 0; j < frameLen; ++j)
		{
			logFilterEnergy[i] += f[j] * melBanks[i][j];
		}
		logFilterEnergy[i] = log10(logFilterEnergy[i]);
	}
	free(melBanks);
	//取对数后作离散余弦变换
	//將上述的 20 個對數能量 Ek帶入離散餘弦轉換，求出 L 階的 Mel- scale Cepstrum 參數，這裡 L 通常取 12。離散餘弦轉換公式如下： 
	//Cm=Sum(cos[m*(k-0.5)*pi/N]*Ek,k=1,n), m=1,2, ..., L 
	//其中 Ek 是由前一個步驟所算出來的三角濾波器和頻譜能量的內積值，N 是三角濾波器的個數。
	for(j = 0;j < DCT_Num; ++j)
	{
		int k;
		mfcc[j] = 0;
		for(k =0;k < TRIANGLE_Num; ++k)
		{
			mfcc[j] += (float)cos(j * (k - 0.5) * PI / TRIANGLE_Num) * logFilterEnergy[k];
		}
	}
	free(logFilterEnergy);
	//加入对数能量
	mfcc[DCT_Num] = logEnergy;
	//标准的倒谱参数MFCC只反映了语音参数的静态特性，语音的动态特性可以用这些静态特征的差分谱来描述。
	//实验证明：把动、静态特征结合起来才能有效提高系统的识别性能
	//差分
	deltaCoeff(mfcc,0, DCT_Num + 1);
	deltaCoeff(mfcc,DCT_Num + 1,2 *( DCT_Num + 1));

	return mfcc;
}
