%% 对于给定采样分析来确定网格的点源执行真空仿真的MATLAB代码

%% 参数初始化
D2 = 0.5;                   % 观测孔径直径
wvl = 1e-6;                 % 光的波长
k = 2 * pi/wvl;              % 波数
Dz = 50e3;                  % 传播距离

DROI = 4*D2;
D1 = wvl*Dz/DROI;
R = Dz;

%% 仿真网格与采样设置
d1 = 10e-3;
d2 = 10e-3;
nscr = 11;
N = 512;
delta1 = d1;                % 源平面网格间距
deltan = d2;                % 观测平面网格间距
n = nscr;                   % 相位屏数量

%% 构造源平面场
% z1=0;
% n_i=1.0;                                % 背景空间折射率
% k0=k*n_i;                               % 背景空间的波数
% ZR=k0*w0^2/2;                           % 瑞利距离
% w = w0 * sqrt(1+(z1/ZR)^2);             % 传播到 z1 处的光束宽度

[x1,y1] = meshgrid((-N/2:N/2-1) * delta1);
[theta,r1] = cart2pol(x1,y1);

pt = exp(-1i * k/(2 * R) * r1.^2)/D1^2 .* sinc(x1/D1) .*sinc(y1/D1) .*exp(-(r1/(4*D1)).^2);

z = (1:n-1) * Dz/(n-1);

sg = exp(-(x1/(0.47 * N * d1)).^16) .* exp(-(y1/(0.47 * N * d1)).^16);
t = repmat(sg,[1,1,n]);
[xn,yn,Uvac] = ang_spec_multi_prop(pt,wvl,delta1,deltan,z,t);
Uvac = Uvac .*exp(-1i * pi /(wvl * R) * (xn.^2+yn.^2));
%% 画图
% 1. 计算辐照度（强度）和相位
Ivac   = abs(Uvac).^2;    % 辐照度分布 W/m^2
Phivac = angle(Uvac);     % 相位分布，单位 rad

% 2. 构造归一化坐标轴（2x/D2, 2y/D2 范围约[-5,5]）
%    注意：xn, yn 单位是米，除以 D2/2 再乘 2 ⇒ 2*xn./D2
x_norm = 2 * xn / D2;
y_norm = 2 * yn / D2;

% (a) 绘制辐照度
figure(1); clf;
imagesc(x_norm(1,:), y_norm(:,1), Ivac);
axis image;                   % 保持 xy 比例
set(gca,'YDir','normal');     % y 轴从下向上
colormap(gray);
cbar = colorbar;
cbar.Label.String = '辐照度 (W/m^2)';
% 根据最大值微调 caxis：示意图中上限约 400，可按需调整
caxis([0, max(Ivac(:))]);
xlabel('2x/D_2');
ylabel('2y/D_2');
title('(a) 辐照度分布');

% (b) 绘制相位
figure(2); clf;
imagesc(x_norm(1,:), y_norm(:,1), Phivac);
axis image;
set(gca,'YDir','normal');
colormap(gray);
cbar = colorbar;
cbar.Label.String = '相位 (rad)';
% 将相位限制在 ±π
caxis([-pi, pi]);
xlabel('2x/D_2');
ylabel('2y/D_2');
title('(b) 相位分布');
