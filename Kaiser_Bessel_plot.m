clc;
clear;

[rt, kf, f] = get_window(51);
figure('NumberTitle', 'off', 'Name', 'Kaiser-Bessel窗口');
title('Kaiser-Bessel窗口');
subplot(1,3,1);
imshow(abs(rt), []);
title('时域');
subplot(1,3,2);
imshow(abs(f), []);
title('频域');
subplot(1,3,3);
imshow(abs(kf),[]);
title('频域的中心3x3区域');

