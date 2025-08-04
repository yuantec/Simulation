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

[pts, TR] = generate_Fibonacci_mesh(N, radius);

% 可视化网格
figure(1);
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
% R = z1(1+(ZR/z1)^2);                         % 等相位面曲率半径
% Phi = atan(z1/ZR);                           % 相位因子
z0 = 0;                                        % 初始面放在 z=0
w0 = 0.02;                                     % 束腰 2 cm

%% 计算初始光场（高斯光束）
U0 = Gaussian_Beam(pts, w0, z0, k);

I0 = abs(U0).^2;
figure(2);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I0, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar;
title('初始高斯光束辐照度分布');

%% 生成初始图拉普拉斯矩阵
L_C = cotangent_Graph_Laplacian(pts, TR);

% 检查对称性与半正定性
assert(issymmetric(L_C), 'Laplacian 矩阵不对称');
eigs_vals = eigs(L_C,1,'smallestabs');
assert(eigs_vals >= -1e-10, 'Laplacian 存在负特征值');
assert(all(abs(sum(L_C,2))<1e-10),'每行和必须为零');

% 特征分解
opts.issym = true;
[Phi, Lambda] = eigs(L_C, N, 'smallestabs', opts);

%% 湍流相关物理量
alpha = 5/3;                  % Kolmogorov 指数
Cn2   = 1e-16;         % 湍流结构常数
r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);   % 平面波 Fried 参量 
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);       % 球面波 Fried 参量
p = linspace(0,Dz,1e3); 
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1))); 

%% 优化求解各相位屏r0scrn
scr_count = 11;       % 相位屏数量 
n = scr_count; 
A = zeros(2,n); 
alpha_vec = (0:n-1)/(n-1); 
A(1,:) = alpha_vec.^(5/3); 
A(2,:) = (1-alpha_vec).^(5/6).*alpha_vec.^(5/6); 
b = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)]; 
x0 = (n/3*r0sw*ones(n,1)).^(-5/3); 
x1 = zeros(n,1); 
rmax = 0.1; 
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:); 
x2(A(2,:)==0) = 50^(-5/3); 
opts = optimoptions('fmincon','Display','none'); 
[X,~,~,~] = fmincon(@(X) sum((A*X - b).^2), x0, [],[],[],[], x1, x2, [], opts); 
r0scrn = X.^(-3/5); 
r0scrn(isinf(r0scrn)) = 1e6; 

%% 生成湍流协方差拉普拉斯矩阵
U = propagate_HalfStep(U0, Phi, Lambda, Dz/2, k); % 初始半步

M_turb = 100;  % 取前 M 模式
nreals = 40;
for idxreal = 1:nreals
    U = U0;                     % 初始场
    for idxscr = 1:n 
        % 湍流相位生成（修正后）
        if idxscr == 1 || r0scrn(idxscr) ~= r0scrn(idxscr-1)
            L_G = construct_Covariance_Laplacian(pts, r0scrn(idxscr), alpha);
            [Psi, Sigma] = eigs(L_G, M_turb, 'largestabs');
        end
        s_n = sqrt(diag(Sigma)/2) .* (randn(M_turb,1) + 1i*randn(M_turb,1));
        S = Psi * s_n;

        % 应用湍流相位
        U = U .* exp(1i * S);

        % 整步传播
        U = propagate_FullStep(U, Phi, Lambda, Dz, k);

        % 自适应网格更新（每2步）
        if mod(idxscr,2) == 0
            [pts, TR] = adaptive_remeshing(pts, U, 30, 1, -30); % θ_max=30°, K=-30dB
            L_C = cotangent_Graph_Laplacian(pts, TR);
            [Phi, Lambda] = eigs(L_C, size(pts,1), 'smallestabs');
        end
    end
end 

%% 可视化
I = abs(U).^2;
figure(3);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I0, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar;
title('湍流光高斯光束辐照度分布');