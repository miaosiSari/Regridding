#include <math.h>     
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "time.h"
#include <sys/time.h>
#define SIZE 500
static void radon(double *pPtr, double *iPtr, double *thetaPtr, int M, int N,    
          int xOrigin, int yOrigin, int numAngles, int rFirst,    
          int rSize);   
   
static char rcs_id[] = "$Revision: 1.13.4.3 $";   
   
#define MAXX(x,y) ((x) > (y) ? (x) : (y))   
   
/* Input Arguments */   
#define I      (prhs[0])   
#define THETA  (prhs[1])   
#define OPTION (prhs[2])   
   
/* Output Arguments */   
#define P      (plhs[0])   
#define R      (plhs[1])   
 
long long gettime(){
	struct timeval s1;
	struct timezone s2;
	gettimeofday(&s1, &s2);
	long long time_microsecond = s1.tv_sec*1000000+s1.tv_usec;
	return time_microsecond;
}

 
void incrementRadon(float *pr, float pixel, float r)   
{   
    int r1;   
    float delta;   
   
    r1 = (int) r;   
    delta = r - r1;   
    pr[r1] += pixel * (1.0 - delta); 
    pr[r1+1] += pixel * delta;  
}   

static void radon(float *pPtr, float *iPtr, float *thetaPtr, int M, int N,    
      int xOrigin, int yOrigin, int numAngles, int rFirst, int rSize)   
{   
    int k, m, n;              /* loop counters */   
    float angle;             /* radian angle value */   
    float cosine, sine;      /* cosine and sine of current angle */   
    float *pr;               /* points inside output array */   
    float *pixelPtr;         /* points inside input array */   
    float pixel;             /* current pixel value */   
    float *ySinTable, *xCosTable;   
    /* tables for x*cos(angle) and y*sin(angle) */   
    float x,y;   
    float r, delta;   
    int r1;   
   
    /* Allocate space for the lookup tables */   
    xCosTable = (float *) calloc(2*N, sizeof(float));  //MATLAB的内存申请函数，对应C语言可以换成calloc函数 
    ySinTable = (float *) calloc(2*M, sizeof(float));   
   
    for (k = 0; k < numAngles; k++) 
	{   
		//每一个theta角，经过radon变化后，就会产生一列数据，这一列数据中，共有rSize个数据
        angle = thetaPtr[k];   
        pr = pPtr + k*rSize;  /* pointer to the top of the output column */   
        cosine = cos(angle);    
        sine = sin(angle);      
   
        /* Radon impulse response locus:  R = X*cos(angle) + Y*sin(angle) */   
        /* Fill the X*cos table and the Y*sin table.                      */   
        /* x- and y-coordinates are offset from pixel locations by 0.25 */   
        /* spaced by intervals of 0.5. */   
		/*
		**radon 变换中，极坐标下，沿r轴的theta角和每一个像素点的分布都是非线性的，而此处采用的是线性radon变换，
		**为了提高精度，把每一个像素点分成其四周四个相邻的像素点来进行计算！x、y坐标的误差是正负0.25
		*/
        for (n = 0; n < N; n++)   
        {   
            x = n - xOrigin;   
            xCosTable[2*n]   = (x - 0.25)*cosine;   //由极坐标的知识知道，相对于变换的原点，这个就是得到了该点的横坐标
            xCosTable[2*n+1] = (x + 0.25)*cosine;   
        }   
        for (m = 0; m < M; m++)   
        {   
            y = yOrigin - m;   
            ySinTable[2*m] = (y - 0.25)*sine;   //同理，相对于变换的原点，得到了纵坐标
            ySinTable[2*m+1] = (y + 0.25)*sine;   
        }   
   
        pixelPtr = iPtr;   
        for (n = 0; n < N; n++)   
        {   
            for (m = 0; m < M; m++)   //便利原矩阵中的每一个像素点
            {   
                pixel = *pixelPtr++;   
                if (pixel != 0.0)   //如果该点像素值不为0，也即图像不连续
                {   
                    pixel *= 0.25;   
					
					//一个像素点分解成四个临近的像素点进行计算，提高精确度
                    r = xCosTable[2*n] + ySinTable[2*m] - rFirst;   
                    incrementRadon(pr, pixel, r);   
   
                    r = xCosTable[2*n+1] + ySinTable[2*m] - rFirst;   
                    incrementRadon(pr, pixel, r);   
   
                    r = xCosTable[2*n] + ySinTable[2*m+1] - rFirst;   
                    incrementRadon(pr, pixel, r);   
   
                    r = xCosTable[2*n+1] + ySinTable[2*m+1] - rFirst;   
                    incrementRadon(pr, pixel, r);   
                }   
            }   
        }   
    }   
                   
    free((void *) xCosTable);   //MATLAB的内存释放函数，对应C语言可以换成free函数
    free((void *) ySinTable);   
}

int main(void){
    int M=SIZE;
	int N=SIZE;
	int xOrigin=((N-1)/2>0)?((N-1)/2):0;
	int yOrigin=((M-1)/2>0)?((M-1)/2):0;
	int temp1=M-1-yOrigin;
	int temp2=N-1-xOrigin;
	int rLast=(int) ceil(sqrt((float) (temp1*temp1+temp2*temp2))) + 1;
	int rFirst=-rLast;
	int rSize=rLast-rFirst+1;
	int numAngles = 181;
	float *thetaPtr= (float*)calloc(numAngles,sizeof(float));
	float *ptr = thetaPtr;
	float deg2rad = 3.141592 / 180.0;
	int k = 0;
	int times = 15;
	for (k = 0; k<numAngles; k++)
		*(ptr++) = numAngles * deg2rad;
    float* I;
    I=(float *)calloc(M*N,sizeof(float));
 
	int p = 0;
	for (p =0;p < M; p++)
	{
		for (k = 0;k < N;k++)
		{
			if (k > M/4 && k < M*3/4)
			{
				if (p > M/4 && p < M*3/4)
					I[p * M + k]=1;
				else
					I[p * M + k]=0;
			}
			else
				I[p * M + k]=0;
		}
	}
	float *cpu_result;
	cpu_result=(float *)calloc(numAngles*rSize,sizeof(float));
	memset(cpu_result,0,numAngles*rSize);
	long long start=0, end=0;
	double sum=0.0, DX=0.0;
    double trial[times+1]={0};
	for (int cnt = 0; cnt < times; ++cnt){
		start = gettime();
		radon(cpu_result, I, thetaPtr, M, N, xOrigin, yOrigin, numAngles, rFirst, rSize);
		end = gettime();
		printf("start=%lld, end=%lld, time=%lld\n", start, end, end-start);
		trial[cnt] = (double)(end-start);
		sum += trial[cnt];
	}
	printf("cpu计算结束\n");
	double average = sum/times;
	printf("Average time=%lf\n", average);
	for (int cnt=0; cnt < times; ++cnt){
        DX += (trial[cnt]-average)*(trial[cnt]-average);
	}		
	printf("DX = %lf\n", DX/times);
	free(I);
	free(thetaPtr);
	free(cpu_result);
	return 0;
}
