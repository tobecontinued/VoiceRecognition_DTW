#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MFCC.h"
//MFCC��Ҫ�ο���Դ
//http://baike.baidu.com/view/2930343.htm#2
//http://blog.csdn.net/c395565746c/article/details/6210920
//http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
//���е������˲����ر���Ҫ
//http://blog.csdn.net/xiaoding133/article/details/8106672?reload

/*
 * �ù��̿��Դﵽ�����򻯽׶ζԾ������ݵ��жϣ�
 * ��Ϊ�������ݵ�ֵ�Ǽ�������ģ�������������Ժ�ֵ���С���ӽ���0������������������ᱣ���ϴ��ֵ
 * ��Ƶ��ĽǶ�����
 * ���������źŵ�Ƶ�ף�ͨ����Ƶ��Խ�߷�ֵԽС���������źŵ�Ƶ����������ʱ���书���׵ķ����½�6dB
 * ��˱���Ը�Ƶ���м��ش���һ���ǽ������ź�ͨ��һ��һ�׸�ͨ�˲���H(z) = 1-0.9375z^-1����ΪԤ�����˲�������Ŀ�����˳���Ƶ����
 * input : sample  ���źŵ�ʱ������
           sampleN ���źŵ�ʱ�����еĳ���
 * output: sample  : Ԥ��Ȩ�غ�������������飬��֤����>=sampleN
 */
void preemphasize(float *sample, int sampleN)   
{     
   int i;   
   float emphFac = (float)0.97;    
   for (i = sampleN - 1; i > 0; --i) 
   {   
        sample[i] = sample[i] - emphFac * sample[i-1];
		//�����ʽ Sr(k) = e(k)(1 - S^-1)���������ӣ���һ�����޼�����Ӧ��ͨ�˲��� H(z) = 1-0.9375z^-1�ǵȼ۵�
		//H(s)  = 1- 0.9375e^-s ����ʽ�任��
		//H(jw) = 1- 0.9375e^-jw 
    }   
    sample[0] = (float)(1.0 - emphFac) * sample[0];
}
/*
 * �������ݵķ�֡
 * �����ź���һ�ֵ��͵ķ�ƽ���źţ����ľ�ֵ����u(x)������غ���R(xl,x2)����ʱ��������ϴ�ı仯��
 * ���о����֣������ź��ڶ�ʱ����Ƶ�����Ա���ƽ�ȣ������ж�ʱƽ�����ԡ�
 * ��ˣ���ʵ�ʴ���ʱ���Խ������źŷֳɺ�С��ʱ���(Լ10~30ms��һ��ȡ֡��20ms��֡��Ϊ֡����1/3~1/2��)����֮Ϊ��֡����
 * ��Ϊ�����źŴ������С��λ��֡��֡�ķ��ص����ֳ�Ϊ֡�ƣ����������źŷֳ�����֡�Ĺ��̳�Ϊ��֡
 * input : sample  :��Ҫ��֡�ĵ�����
		   sampleN ���źŵ�ʱ�����еĳ���
		   frameLen ��ÿһ֡�ĳ��ȣ����Ϊ2���ݣ�����FFT,DCT�ȼ���
		   frameOffest : frameLen * 0.33
   output : pframeN ��֡������������
   return : �ú�������һ������[֡��][֡��]����ʾ��֡�������,
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
 * ������������Ϊ�˱��������źŵĶ�ʱƽ���ԣ����ô������������ɽضϴ����µ�GibbsЧӦ
 * ����������:
      Whm(n) = 0.5-0.46*cos(2*PI*n/(N-1)) ((0��n��N-1))
	      = 0                             (Other)
 * GibbsЧӦ:�����в�������ĺ���,ȡ����Ҷ�����󣬼�ʹ��ȡ������������ڲ��������Բ�����
			 ��ԭ��������Խ��㸽���Ĳ��Σ����ǲ��ɱ���Ĵ���������񵴣��Ӷ���Խ��㸽��ĳЩ
			 ��ĺ���ֵ����1,���γɴӹ��壬���ż�����ȡ���������࣬����񵴴��ڵ�ʱ�佫���̣�
			 ������ֵ����Խ9%�Ĺ̶�ֵ--���ź�������ϵͳ����P104
 * ������������ԭ����ʱ�����
 * input  : frame :��Ҫ�Ӵ���һ֡����
          ��windowLen ���ĳ��ȣ�ͨ��ȡ֡��
 * output : frame �Ӵ��������
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
 * ���ٸ���Ҷ�任 --���ź�������ϵͳ���ھ���
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
* ��fsΪʱ��ĳ���Ƶ��
*  FFT��Fo��Ƶ��ĳ������
* �� Fo = fs/n	
	
 * input :  Real :ʱ�����е�ʵ��
            Imag :ʱ�����е��鲿
		    n    :ʱ�����еĳ���,nӦ������ n == 2^k;
			������ʱ��0ֵ�ﵽ���г���Ҫ��       
 * output : Real :����Ҷ�任���ʵ��
            Imag :����Ҷ�任����鲿

* Wnk = e^(-j*2*pi*k/n)
*/
//����x�洢ʱ�����е�ʵ��������y�洢ʱ�����е��鲿
//n����N��FFT��sign=1ΪFFT��sign=-1ΪIFFT
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
	//�������ͼ�������±꣨��λ������
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
	//����ÿһ���������lΪĳһ����iΪͬһ���Ĳ�ͬȺ��ʹ��ͬһ�ڴ棨��λ���㣩
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
	//�������IFFT���ٳ���N
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
 * ���Ǵ�ͨ�˲�������Ƶ���������
 * ��o,c,h�ֱ�Ϊ�����˲��������ޣ����ģ�����Ƶ��y
 * on��cn��hnΪƵ�ʳ���Ƶ��������
 * ���� W(k) = (k - on��/(cn-on) on < = k  < =cn
			   (hn-k) /(hn- cn)  cn < = k  < =hn	
 * input  : F	 :��Ҫ�˲�����
            n	 :�׵ĳ���
			Fbeta   : Ƶ������2���Ƶ�ʲ�
		    oFreq :����Ƶ��
		    hFreq :����Ƶ��
 * output : v 
 */
void TriangleFilter(float *F,int Fbeta,int n,float oFreq,float hFreq)
{
	/*
	 * ��fsΪʱ��ĳ���Ƶ��
	 *  FFT��Fo��Ƶ�ڵĳ������
	 * �� Fo = fs/n
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
 * �����˲������ǵķ�Χ���������˶���һ���ٽ����,�Դ��� ģ���˶����ڱ�ЧӦ(������˵Ӧ����Ƶ����ڱ�ЧӦ)
 * �Թ������˲�
 * �����Ƶ��Ϊfs���ɳ��������źŵ����Ƶ��ΪfgΪfs/2(��fg = fs/n *(n-1))
 * ����������fmel = 2595 *log(1 + f /700hz)
 * fΪƵ��
 * �õ���ߵ��������� fmax
 * ��Mel�̶ȷ�Χ�ڣ����������˲���������Ƶ������ȼ�������Էֲ�
 * ���Լ����������������˲���������Ƶ�ʵļ��Ϊ
 * deta = fmax / (k+1)
 * ÿһ���������˲���������Ƶ��c(l) __��MelƵ�����ϵȼ���ֲ�__����o(l),c(l),h(l) �ֱ��ǵ�l ���������˲��������ޣ����ģ�������Ƶ�ʣ�
 * �������������˲���֮������ޣ����ģ�����Ƶ�ʵĹ�ϵ���£�c(l)=h(l-1)=o(l+1)
 * input : F :��Ҫ�˲��Ĺ�����
           n :����
		   SampleRate��������
		   k :�˲�������
 * return : ���س���ΪK�ĵ����飬ÿ��ֵΪͨ��ÿ���˲����Ķ������������x��Ӎ̖��ƽ���ͣ���ȡ�� 10 ��׵Č���ֵ���ٳ��� 10
 * ���Ӧ��Ϊ k������
 * ��K����Ϊͨ����K���˲������������
 * Ŀǰ�˺�������
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
 * �����˲�����
 * input: k : �˲�����ĸ���
 *		  n : ��Ҫ�˲������г���
 *	sampleRate : ������
 * return �� �˲����������У�Ϊfloat[k][n]����
*/
float** melbank(int k, int n, int sampleRate)
{
	int i,j;
	float melDelta = freqtoMel(sampleRate / 2.0) / (k + 1);//�˲���������Ƶ�ʼ��
	float freqDelta = sampleRate / n;
	float freq;
	float *f = (float*)calloc(k + 2,sizeof(float));//�˲����������Ƶ��
	float **melbanks = (float**)calloc(k,sizeof(float*));
		
	for (i = 0; i < k; ++i)
	{
		melbanks[i] = calloc(n,sizeof(float));
		for (j = 0; j < n; ++j)
			melbanks[i][j] = 0.0;
	}
	//���˲�����Ƶ�������Ƶ�ʣ�����mel��Ϊ�Ⱦ�ֲ�
	f[0] = 0.0;
	for (i = 1; i <= k + 1; ++i)
		f[i] = f[i - 1] + melDelta;
	for (i = 1; i <= k + 1; ++i)
		f[i] = meltoFreq(f[i]);
	//�˲�����
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
 * ��m[beg,end)�Ĳ��,�������
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
 * ÿ֡�������ټ�12;Ϊdct��������13άMFCC
 * �ڼ���һ�β�֣��Ͷ��β�֣�����39ά
 */
float* FrametoMFCC(float* f,int frameLen,int SampleRate)
{
	int i,j;
	float *fImag,*mfcc,
		  logEnergy = 0.0,
		  *logFilterEnergy;
	float **melBanks;
	mfcc = (float*)calloc(MFCC_DIMENSION,sizeof(float));
	//�Ӵ�
	HammingWindowing(f ,frameLen);
	//����һ֡�Ķ������� ���x��һ�������Ӎ̖��ƽ���ͣ���ȡ�� 10 ��׵Č���ֵ���ٳ��� 10
	for(j = 0;j < frameLen;++j)
	{
		logEnergy += f[j] * f[j];
	}

	logEnergy = (float)log10(logEnergy) * 10;
	//���㹦����
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
	//��������ͨ���˲�����,����ÿ���˲����Ķ������� 
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
	//ȡ����������ɢ���ұ任
	//�������� 20 ���������� Ek�����xɢ�N���D�Q����� L �A�� Mel- scale Cepstrum �������@�e L ͨ��ȡ 12���xɢ�N���D�Q��ʽ���£� 
	//Cm=Sum(cos[m*(k-0.5)*pi/N]*Ek,k=1,n), m=1,2, ..., L 
	//���� Ek ����ǰһ�����E�����������ǞV�������l�V�����ăȷeֵ��N �����ǞV�����Ă�����
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
	//�����������
	mfcc[DCT_Num] = logEnergy;
	//��׼�ĵ��ײ���MFCCֻ��ӳ�����������ľ�̬���ԣ������Ķ�̬���Կ�������Щ��̬�����Ĳ������������
	//ʵ��֤�����Ѷ�����̬�����������������Ч���ϵͳ��ʶ������
	//���
	deltaCoeff(mfcc,0, DCT_Num + 1);
	deltaCoeff(mfcc,DCT_Num + 1,2 *( DCT_Num + 1));

	return mfcc;
}
