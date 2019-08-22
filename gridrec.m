clc;
clear;
%��������ָ�����Ķ���:
%Regridding reconstruction algorithm for real-time
%tomographic imaging
N = 331;%Ŀǰ���NΪ����
dtheta = 0.1;%ÿ��0.1��ͶӰ
thetas = 0:dtheta:180-dtheta;
len = length(thetas);
p = phantom(N);%p��һ��NxN��Shepp-Logan Phantom
%%ͶӰ��Zero-Padding�ļ���
if mod(N, 2) == 1
    padding = N - 1;
else
    padding = N;
end
isshow = 1;%�Ƿ�չʾKaiser-Besselͼ=
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
%Radon�任
for cnt = 1:len
    theta = thetas(cnt);
    p_rot = imrotate(p, -theta, 'bilinear', 'crop');%����Bicubic����ת��ֵ���ܴ�������Ϊ������Ŀ����
    r1(:, cnt) = sum(p_rot, 2);%���з���Ϊx�ᣬ�з���Ϊy���y��͡�
end
for cnt = 1:N
    %���0��ʹ�ò���Ƶ�ʴ��ڵ����ο�˹��Ƶ�ʣ��������������2.3.1
    r(cnt+halfpadding, :) = r1(cnt, :);
end


R = fftshift(fft(ifftshift(r)));
Q = zeros(N+padding);

%���ڸĽ���ramp filter, �������2.3.2 Constant Offset
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
for cnt = 1:N+padding %��wѭ��
    w = cnt-(N+padding+1)/2; 
    for t = 1:len % ��thetaѭ��
       theta = thetas(t);
       theta_arc = theta*pi/180;
       costheta = cos(theta_arc);
       sintheta = sin(theta_arc);
       P_w_theta = R(cnt, t);%Fourier Slice,��һ���൱������ʽ(6)�е�P��
       wcostheta = w*costheta;
       wsintheta = w*sintheta;
       %ֻ��Ҫ�Ծ���˵Ĵ�С����ѭ��
       %���ھ������Kaiser-Bessel���ھ���fftshift(fft2(ifftshift()))�õ���������ԭ��Ҳ���м�
       for height = -halfkh:1:halfkh
           for width = -halfkw:1:halfkw
               Cx = round(wcostheta+(N+padding+1)/2+height);
               Cy = round(wsintheta+(N+padding+1)/2+width);
               %Cx, Cy�൱��ʽ(6)�е�U, V�����ھ���˿�Ⱥ�С���������ѭ���е�ÿ��w, \theta, ֻ��Ҫ��
               %w cos\theta, w sin\theta�����ĺ;����һ���������ѭ������
               %ע�⵽����Ҷ�任�������ԣ���Ҫ��һЩ�߽��������д���
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
               Cvalue = kf(height+halfkh+1, width+halfkw+1); %�൱��ʽ(6)�е�C
               %fprintf('Cx=%d, Cy=%d\n',Cx,Cy);
               Q(Cx, Cy) = Q(Cx, Cy)+ Cvalue * P_w_theta * RAMP(cnt);
           end
       end
    end
end
Q = Q*1/N*dtheta/N;
q = fftshift(ifft2(ifftshift(Q)));%��任
p_recon = real(q./rt);%����Kaiser-Bessel����
p_recon(p_recon>1) = 1.0;
p_recon(p_recon<0) = 0.0;
figure;
p_recon = p_recon(1+halfpadding:N+halfpadding, 1+halfpadding:N+halfpadding);%����padding
toc;
imshow(p_recon,[]);
title('Gridrec');