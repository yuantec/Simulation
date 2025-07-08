clc; clear; close all;
%% 基本物理与系统参数 
D2 = 0.5;            % 观测孔径直径 (m)
wvl = 1e-6;          % 波长 (m)
k = 2*pi/wvl;        % 波数
Dz = 50e3;           % 总传播距离 (m)

%% 初始网格与采样设置
N = 512;                            % 顶点数量
d1 = 10e-3;                         % 源面采样间距 (m)
d2 = 10e-3;                         % 观测面采样间距 (m)
radius = (N*d1)/sqrt(pi);           % 圆盘半径 ≈ (N d1)/√π
nscr = 11;
n = nscr;                           % 相位屏数量
[pts, TR] = generate_Fibonacci_mesh(N, radius);

% 可视化网格
figure;
triplot(TR);
axis equal;
title('初始 Fibonacci Delaunay 网格');
xlabel('x (m)');
ylabel('y (m)');

%% 普通高斯光束参数
% z1 = 0;
% n_i = 1.0;                                   % 背景空间折射率
% k0 = k*n_i;                                  % 背景空间的波数
% ZR = k0*w0^2/2;                              % 瑞利距离                                        
% w = w0*sqrt(1+(z1/ZR)^2);                    % 传播到z处的束宽
% R = z1(1+(ZR/z1)^2);                        % 等相位面曲率半径
% Phi = atan(z1/ZR);                           % 相位因子
z0 = 0;                  % 初始面放在 z=0
w0 = 0.02;               % 束腰 2 cm

%% 计算初始光场（高斯光束）
U0 = Gaussian_Beam(pts, w0, z0, k);

I0 = abs(U0).^2;
figure;
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I0, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar;
title('初始高斯光束辐照度分布');
%% 生成初始图拉普拉斯矩阵
L_C = cotangent_Graph_Laplacian(pts, TR);

% 检查对称性与半正定性
assert(issymmetric(L_C), 'Laplacian 矩阵不对称');
eigs_vals = eigs(L_C,1,'smallestabs');
assert(eigs_vals >= -1e-10, 'Laplacian 存在负特征值');

% 特征分解
opts.issym = true;
[Phi, Lambda] = eigs(L_C, N, 'smallestabs', opts);

%% 真空传播
dz = Dz / n;  % 总传播距离分成n段
for step = 1:n
    U0 = propagate_HalfStep(U0, Phi, Lambda, dz, k);
    
    U0 = propagate_HalfStep(U0, Phi, Lambda, dz, k);
end

% 可视化结果
figure;
scatter(pts(:,1), pts(:,2), 15, abs(U0).^2,'filled');
axis equal; colorbar;
title('真空传播后辐照度分布');
