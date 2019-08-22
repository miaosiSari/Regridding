clc;
clear all;

N = 331;
dtheta = 0.1;%每隔0.1度投影
thetas = 0:dtheta:180-dtheta;
len = length(thetas);
p = phantom(N);%p是一个NxN的Shepp-Logan Phantom

r = radon(p, thetas);
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
R = fftshift(fft2(ifftshift(rfill)));
Q = zeros(N1+padding);

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
for cnt = 1:N1+padding %对w循环
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
               Cvalue = kf(height+halfkh+1, width+halfkw+1); %相当于式(6)中的C
               %fprintf('Cx=%d, Cy=%d\n',Cx,Cy);
               Q(Cx, Cy) = Q(Cx, Cy)+ Cvalue * P_w_theta * RAMP(cnt);
           end
       end
    end
end
Q = Q*1/N1*dtheta/N1;
q = fftshift(ifft2(ifftshift(Q)));%逆变换
p_recon = real(q./rt);%消除Kaiser-Bessel窗口
p_recon(p_recon>1) = 1.0;
p_recon(p_recon<0) = 0.0;
%p_recon = p_recon(1+halfpadding:N+halfpadding, 1+halfpadding:N+halfpadding);%消除padding
toc;
figure;
imshow(p_recon,[]);
title('Gridrec');