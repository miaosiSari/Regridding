function [raw_kaiser_bessel, f_cropped, f] = get_window(N, crop)
%%
%N: Shepp-Logan Phantom的尺寸
%取Kaiser-Bessel窗口的中心(2*crop+1, 2*crop+1)的能量集中部分
%根据论文Regridding reconstruction algorithm for real-time
%tomographic imaging公式(6), 计算每个频域点Q(U,
%V)需要NM次计算(M为投影数量),而一共有N^2个频域点，因此需要N^3M次计算
%Gridrec的巧妙之处在于使用了一种Kaiser-Bessel窗口,它在时域上能量分散且处处不为0(这样就可以将Kaiser窗口消去),而在频域上能量集中在直流分量(DC)附近,大部分区域为0
%如此,只需要循环直流分量附近3x3的非零区域即可!
if ~exist('crop', 'var')
    crop = 1; 
end
raw_kaiser_bessel = w_matrix(N);%原始的Kaiser-Bessel窗口K
%w_matrix:生成NxN的Kaiser-Bessel窗口
%w_matrix调用w_vector, w_vector生成长度为n的Kaiser-Bessel向量
%w_vector调用w_single, 生成单个Kaiser-Bessel采样点，计算公式位于:
f = fftshift(fft2(ifftshift(raw_kaiser_bessel)));%K的傅里叶变换
if mod(N, 2) == 1
    f_cropped = f((N+1)/2-crop:(N+1)/2+crop,(N+1)/2-crop:(N+1)/2+crop);
else
    f_cropped = f(N/2-crop:(N+1)/2+crop, (N+1)/2-crop:(N+1)/2+crop);%截取
end