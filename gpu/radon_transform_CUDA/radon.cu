#include <stdio.h>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "time.h"
#include <sys/time.h>
#define SIZE 500   
__global__ static void radon_cuda_core(float *gpuOutput,float *gpuInput,float *gpuAngles,int M,int N,int xOrgin,int yOrgin,int numAngles,int rFirst,int rSize)
{
	const int tid=threadIdx.x;	
	float angle=gpuAngles[tid];
	float *pOutput=gpuOutput+tid*rSize;		
	float sine=sin(angle);
	float cosine=cos(angle);
	int m,n;
 
	for(m=0;m<rSize;m++)
		pOutput[m]=0.0;
 
	float *pInput=gpuInput;
	for(n=0;n<N;n++)
	{
		for(m=0;m<M;m++)
		{
			float value=*pInput++;
			if(value!=0.0)
			{
				value*=0.25;
				float x=n-xOrgin;
				float y=yOrgin-m;
 
				int r1;
				float delta;
					
				float r=(x-0.25)*cosine+(y-0.25)*sine-rFirst;
				r1=(int)r;
				delta=r-r1;
				pOutput[r1]+=value*(1.0-delta);
				pOutput[r1+1]+=value*delta;
 
				r=(x-0.25)*cosine+(y+0.25)*sine-rFirst;
				r1=(int)r;
				delta=r-r1;
				pOutput[r1]+=value*(1.0-delta);
				pOutput[r1+1]+=value*delta;
 
				r=(x+0.25)*cosine+(y+0.25)*sine-rFirst;
				r1=(int)r;
				delta=r-r1;
				pOutput[r1]+=value*(1.0-delta);
				pOutput[r1+1]+=value*delta;
 
				r=(x+0.25)*cosine+(y-0.25)*sine-rFirst;
				r1=(int)r;
				delta=r-r1;
				pOutput[r1]+=value*(1.0-delta);
				pOutput[r1+1]+=value*delta;
			}
		}
	}	
}
static void radon_cuda(float *pPtr, float *iPtr, float *thetaPtr, int M, int N,    
	int xOrigin, int yOrigin, int numAngles, int rFirst, int rSize)
{
	float *gpuInput;
	float *gpuOutput;
	float *gpuAngles;
	cudaMalloc((void **)&gpuInput,sizeof(float)*M*N);
	cudaMalloc((void **)&gpuOutput,sizeof(float)*numAngles*rSize);
	cudaMalloc((void **)&gpuAngles,sizeof(float)*numAngles);
	cudaMemcpy(gpuInput,iPtr,sizeof(float)*M*N,cudaMemcpyHostToDevice);
	cudaMemset(gpuOutput,0,numAngles*rSize);
	cudaMemcpy(gpuAngles,thetaPtr,sizeof(float)*numAngles,cudaMemcpyHostToDevice);
 
	radon_cuda_core<<<1,numAngles,0>>>(gpuOutput,gpuInput,gpuAngles,M,N,xOrigin,yOrigin,numAngles,rFirst,rSize);
 
	cudaMemcpy(pPtr,gpuOutput,sizeof(float)*numAngles*rSize,cudaMemcpyDeviceToHost);
 
	cudaFree(gpuInput);
	cudaFree(gpuOutput);
	cudaFree(gpuAngles);
}

long long gettime(){
	struct timeval s1;
	struct timezone s2;
	gettimeofday(&s1, &s2);
	long long time_microsecond = s1.tv_sec*1000000+s1.tv_usec;
	return time_microsecond;
}
 
int main()
{
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
	float* I=(float *)calloc(M*N,sizeof(float));
 
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
	float *gpu_result;
	gpu_result=(float *)calloc(numAngles*rSize,sizeof(float));
	memset(gpu_result,0,numAngles*rSize);
	long long start=0, end=0;
	double sum=0.0, DX=0.0;
    double trial[times+1]={0};
	for (int cnt = 0; cnt < times; ++cnt){
		start = gettime();
		radon_cuda(gpu_result, I, thetaPtr, M, N, xOrigin, yOrigin, numAngles, rFirst, rSize);
		end = gettime();
		printf("start=%lld, end=%lld, time=%lld\n", start, end, end-start);
		trial[cnt] = (double)(end-start);
		sum += trial[cnt];
	}
	double average = sum/times;
	printf("Average time=%lf\n", average);
	for (int cnt=0; cnt < times; ++cnt){
        DX += (trial[cnt]-average)*(trial[cnt]-average);
        printf("%lld ",(long long)trial[cnt]);
	}		
        printf("\n");
	printf("DX = %lf\n", DX/times);
	free(I);
	free(thetaPtr);
	free(gpu_result);
	return 0;
}
