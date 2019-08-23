clc;
clear;

[rt, kf, f] = get_window(51);
figure('NumberTitle', 'off', 'Name', 'Kaiser-Bessel����');
title('Kaiser-Bessel����');
subplot(1,3,1);
imshow(abs(rt), []);
title('ʱ��');
subplot(1,3,2);
imshow(abs(f), []);
title('Ƶ��');
subplot(1,3,3);
imshow(abs(kf),[]);
title('Ƶ�������3x3����');

