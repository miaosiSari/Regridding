clc;
clear all;

N = 331;
dtheta = 0.1;%ÿ��0.1��ͶӰ
thetas = 0:dtheta:180-dtheta;
len = length(thetas);
p = phantom(N);%p��һ��NxN��Shepp-Logan Phantom

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

%���ڸĽ���ramp filter, �������2.3.2 Constant Offset
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
for cnt = 1:N1+padding %��wѭ��
    w = cnt-(N1+padding+1)/2; 
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
               Cx = round(wcostheta+(N1+padding+1)/2+height);
               Cy = round(wsintheta+(N1+padding+1)/2+width);
               %Cx, Cy�൱��ʽ(6)�е�U, V�����ھ���˿�Ⱥ�С���������ѭ���е�ÿ��w, \theta, ֻ��Ҫ��
               %w cos\theta, w sin\theta�����ĺ;����һ���������ѭ������
               %ע�⵽����Ҷ�任�������ԣ���Ҫ��һЩ�߽��������д���
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
               Cvalue = kf(height+halfkh+1, width+halfkw+1); %�൱��ʽ(6)�е�C
               %fprintf('Cx=%d, Cy=%d\n',Cx,Cy);
               Q(Cx, Cy) = Q(Cx, Cy)+ Cvalue * P_w_theta * RAMP(cnt);
           end
       end
    end
end
Q = Q*1/N1*dtheta/N1;
q = fftshift(ifft2(ifftshift(Q)));%��任
p_recon = real(q./rt);%����Kaiser-Bessel����
p_recon(p_recon>1) = 1.0;
p_recon(p_recon<0) = 0.0;
%p_recon = p_recon(1+halfpadding:N+halfpadding, 1+halfpadding:N+halfpadding);%����padding
toc;
figure;
imshow(p_recon,[]);
title('Gridrec');