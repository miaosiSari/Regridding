clc;
clear all;

N = 521;
dtheta = 0.1;%每隔0.1度投影
thetas = 0:dtheta:180-dtheta;
len = length(thetas);
trueimage = 0;
if trueimage == 0 %不使用真实图像而使用Shepp-Logan Phantom
    p = phantom(N);%p是一个NxN的Shepp-Logan Phantom
else
    p = im2double(rgb2gray(imread('1.png')));
    p = imresize(p, [N, N]);
end

%需要注意radon函数以行为x轴正向,列为y轴正向, 坐标系为左手坐标。
%因此最后得到的图像需要用flipud函数(flip up and down)作颠倒。
%和gridrec.m不同,gridrec_radon使用radon函数作投影。如果使用imrotate是N根射线，如果使用radon是
% 2*ceil(norm(size(p)-floor((size(p)-1)/2)-1))+3根射线(请打开radon的定义)
%事实证明投影射线越多效果越好
[r, xp] = radon(p, thetas);
N1 = size(r, 1);
padding = N1-1;
halfpadding = padding/2;
rfill = zeros(N1+padding, len);
[rt, kf, f] = get_window(N1+padding);
[kh, kw] = size(kf);
halfkh = (kh-1)/2;
halfkw = (kw-1)/2;
for cnt = 1:N1
    rfill(cnt+halfpadding, :) = r(cnt, :);
end
R = fftshift(fft(ifftshift(rfill)));
Q = zeros(N1+padding);%在Matlab中zeros(k)实际上创建的是零填充的二维数组NxN
c = N1+padding-N;

%关于改进的ramp filter, 详见论文2.3.2 Constant Offset
ramp = zeros(N1+padding, 1);
for cnt = 1:N1+padding
    realT = cnt-(N1+padding+1)/2;
    if realT == 0
       ramp(cnt) = 1/4; 
    elseif mod(realT, 2) == 1
       ramp(cnt) = -1/(realT*pi)^2;
    else
       ramp(cnt) = 0;
    end
end
RAMP = fftshift(fft(ifftshift(ramp)));
tic;
for cnt = 1 : N1+padding %对w循环
    w = cnt-(N1+padding+1)/2; 
    for t = 1:len % 对theta循环
       theta = thetas(t);
       theta_arc = theta*pi/180;
       costheta = cos(theta_arc);
       sintheta = sin(theta_arc);
       P_w_theta = R(cnt, t);%Fourier Slice,这一项相当于论文式(6)中的P。
       wcostheta = w*costheta;
       wsintheta = w*sintheta;
       %只需要对卷积核的大小进行循环
       %由于卷积核是Kaiser-Bessel窗口经过fftshift(fft2(ifftshift()))得到，其坐标原点也在中间
       for height = -halfkh:1:halfkh
           for width = -halfkw:1:halfkw
               Cx = round(wcostheta+(N1+padding+1)/2+height);
               Cy = round(wsintheta+(N1+padding+1)/2+width);
               %Cx, Cy相当于式(6)中的U, V。由于卷积核宽度很小，对于外层循环中的每个w, \theta, 只需要对
               %w cos\theta, w sin\theta附近的和卷积核一样大的邻域循环即可
               %注意到傅里叶变换的周期性，需要对一些边界条件进行处理
               if Cx > N1 + padding
                   Cx = Cx - (N1 + padding);
               end
               if Cx < 1
                   Cx = Cx + (N1 + padding);
               end
               if Cy > N1 + padding
                   Cy = Cy - (N1 + padding);
               end
               if Cy < 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                   Cy = Cy + (N1 + padding);
               end
               Cvalue = kf(width+halfkw+1, height+halfkh+1); %相当于式(6)中的C
               %fprintf('Cx=%d, Cy=%d\n',Cx,Cy);
               Q(Cy, Cx) = Q(Cy, Cx)+ Cvalue * P_w_theta * RAMP(cnt);
           end
       end
    end
end
Q = Q*1/(N1+padding)*dtheta/(N1+padding);%这里系数(N1+padding)其实并非很重要，因为图像信息由香味决定
q = fftshift(ifft2(ifftshift(Q)));%逆变换
p_recon = real(q./rt);%消除Kaiser-Bessel窗口
p_recon(p_recon>1) = 1.0;
p_recon(p_recon<1e-3) = 0.0;
disp(mean(p_recon(:)));
disp(mean(p(:)));
if mod(c, 2) == 0
   p_recon = p_recon(1+c/2:end-c/2, 1+c/2:end-c/2);%消除padding
else
   p_recon = p_recon(1+floor(c/2):end-ceil(c/2), 1+floor(c/2):end-ceil(c/2));
end%需要截取中间区域。如果不截取，图像边缘会有亮斑，灰度值接近1。在伸缩变换下有效区域会接近0
p_recon = flipud(p_recon);
toc;
figure;
subplot(1,2,1);
imshow(p,[]);
title('Original');
subplot(1,2,2);
imshow(p_recon,[]);
title('Gridrec Reconstruction');
p_recon = p_recon*mean(p(:))/mean(p_recon(:));
imwrite(p_recon, 'recon.png');
imwrite(p, 'p.png');
[psnr, ssim] = calc_psnr_ssim('p.png', 'recon.png'); %PSNR=32.8999, SSIM=0.9649, excellent!
%但是低分辨率重建很模糊，怎么办?
%难道Deblur?