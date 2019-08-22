clc;
clear;
%代码中所指的论文都是:
%Regridding reconstruction algorithm for real-time
%tomographic imaging
N = 331;%目前针对N为奇数
dtheta = 0.1;%每隔0.1度投影
thetas = 0:dtheta:180-dtheta;
len = length(thetas);
p = phantom(N);%p是一个NxN的Shepp-Logan Phantom
%%投影的Zero-Padding的计算
if mod(N, 2) == 1
    padding = N - 1;
else
    padding = N;
end
isshow = 1;%是否展示Kaiser-Bessel图=
halfpadding = padding/2;
r1 = zeros(N, len);
r = zeros(N+padding, len);
[rt, kf, f] = get_window(N+padding);%rt: Raw Time Domain, kf: Cropped Frequency Domain, f:fft2 of rt.
if isshow
   figure;
   subplot(1,3,1);
   imshow(abs(rt), []);
   title('Kaiser-Bessel Window.');
   subplot(1,3,2);
   imshow(abs(kf), []);
   title('Cropped Kaiser-Bessel Window (Frequency Domain)');
   subplot(1,3,3);
   imshow(abs(f), []);
   title('Original Kaiser-Bessel Window (Frequency Domain), energy gathers at the DC');
end
[kh, kw] = size(kf);
halfkh = (kh-1)/2;
halfkw = (kw-1)/2;
%Radon变换
for cnt = 1:len
    theta = thetas(cnt);
    p_rot = imrotate(p, -theta, 'bilinear', 'crop');%这种Bicubic的旋转插值可能带来误差，因为射线数目不足
    r1(:, cnt) = sum(p_rot, 2);%以行方向为x轴，列方向为y轴对y求和。
end
for cnt = 1:N
    %填充0，使得采样频率大于等于奈奎斯特频率，抗混叠。见论文2.3.1
    r(cnt+halfpadding, :) = r1(cnt, :);
end


R = fftshift(fft(ifftshift(r)));
Q = zeros(N+padding);

%关于改进的ramp filter, 详见论文2.3.2 Constant Offset
ramp = zeros(N+padding, 1);
for cnt = 1:N+padding
    realT = cnt-(N+padding+1)/2;
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
for cnt = 1:N+padding %对w循环
    w = cnt-(N+padding+1)/2; 
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
               Cx = round(wcostheta+(N+padding+1)/2+height);
               Cy = round(wsintheta+(N+padding+1)/2+width);
               %Cx, Cy相当于式(6)中的U, V。由于卷积核宽度很小，对于外层循环中的每个w, \theta, 只需要对
               %w cos\theta, w sin\theta附近的和卷积核一样大的邻域循环即可
               %注意到傅里叶变换的周期性，需要对一些边界条件进行处理
               if Cx > N + padding
                   Cx = Cx - (N + padding);
               end
               if Cx < 1
                   Cx = Cx + (N + padding);
               end
               if Cy > N + padding
                   Cy = Cy - (N + padding);
               end
               if Cy < 1
                   Cy = Cy + (N + padding);
               end
               Cvalue = kf(height+halfkh+1, width+halfkw+1); %相当于式(6)中的C
               %fprintf('Cx=%d, Cy=%d\n',Cx,Cy);
               Q(Cx, Cy) = Q(Cx, Cy)+ Cvalue * P_w_theta * RAMP(cnt);
           end
       end
    end
end
Q = Q*1/N*dtheta/N;
q = fftshift(ifft2(ifftshift(Q)));%逆变换
p_recon = real(q./rt);%消除Kaiser-Bessel窗口
p_recon(p_recon>1) = 1.0;
p_recon(p_recon<0) = 0.0;
figure;
p_recon = p_recon(1+halfpadding:N+halfpadding, 1+halfpadding:N+halfpadding);%消除padding
toc;
imshow(p_recon,[]);
title('Gridrec');