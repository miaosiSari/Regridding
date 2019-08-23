clear all;
clc;

%======================================================================
%In this program we verify the case that the size of Shepp-Logan is odd
%The even case is similar.
l = 129;
M = 2*l+1;
%======================================================================
x = phantom(M);
theta = 20;
X = fftshift(fft2(ifftshift(x)));
y = imrotate(x, -theta, 'bilinear', 'crop');
z = sum(y, 2);
Z = fftshift(fft(ifftshift(z)));
W = zeros(size(Z));
theta_arc = theta*pi/180;
costheta = cos(theta_arc);
sintheta = sin(theta_arc);
for cnt = 1:M
    omega = cnt-(M+1)/2;
    x = omega*costheta+(M+1)/2;
    y = omega*sintheta+(M+1)/2;
    newupx = floor(x);
    newdownx = ceil(x);
    newlefty = floor(y);
    newrighty = ceil(y);
    W(cnt) = X(newupx, newlefty)*(1-x+newupx)*(1-y+newlefty)+...
        +X(newupx, newrighty)*(1-x+newupx)*(y-newlefty)+...
        +X(newdownx, newlefty)*(x-newupx)*(1-y+newlefty)+...
        +X(newdownx, newrighty)*(x-newupx)*(y-newlefty);
    fprintf('[M=%f] omega=%f, x=%f, y=%f,\n newupx=%f, newdownx=%f, newlefty=%f, newrighty=%f\n\n',M, omega, x, y, newupx, newdownx, newlefty, newrighty);
end
z1 = fftshift(ifft2(ifftshift(Z)));
w1 = fftshift(ifft2(ifftshift(W)));
w2 = real(w1);
error = sum((w2-z).^2);
signal = sum(z.^2);
SNR = 10*log10(signal/error);
plot(z);
grid on;hold on; 
plot(w2);
legend('FFT of projections','FFT2 of the original image');