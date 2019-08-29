#radon_transform: A CUDA implementation of the Radon transform.

##How to run:
For CPU version: ./radon.cpu
For GPU version ./radon.gpu

##How to compile:
g++ -o radon.cpu radon.cpp (CPU Version)
nvcc -o radon.gpu radon.cu (GPU Version)

The source code owes to Matlab (CPU Version) and https://blog.csdn.net/celte/article/details/9826505 (CUDA Version).
