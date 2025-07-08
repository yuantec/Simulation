D = 2;               % 相位屏长度
r0 = 0.1;           % 大气相干长度
N = 256;            % 采样点数量
L0 = 100;           % 外尺度
l0 = 0.01;          % 内尺度

delta = D/N;        % 采样间隔
X = (-N/2:N/2-1) * delta;
Y = X;
[phz_lo,phz_hi] = ft_sh_phase_screen(r0,N,delta,L0,l0);
phz = phz_lo + phz_hi;


%%% 二维图像 %%%
figure (5)
pcolor(X,Y,phz); 
colormap("jet")
axis on;
shading interp;
axis square;  %输出正方形图像
colorbar;
title("低频补偿后的大气湍流随机相位屏（二维）");
% xlabel('{\itx}(m)');
% ylabel('{\ity}(m)');
%%% 二维图像 %%%
%%% 三维图像 %%%
figure (55)
surf(X,Y,phz); 
colormap("jet")
axis on;
shading interp;
axis square;  %输出正方形图像
colorbar;
title("低频补偿后的大气湍流随机相位屏(三维)");
%%% 三维图像 %%%