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
% 预先定义
M_turb = 100;                            % KL 模式数量
dz = Dz / scr_count;                     % 每段距离
opts.issym = true;

% 初始半步衍射
U = propagate_HalfStep(U0, Phi, Lambda, dz/2, k);

for idxreal = 1:nreals
    % 重置场为初始半步后状态
    U_current = U;

    for idxscr = 1:scr_count
        %% —— 湍流相位屏 （KL展） —— 
        % 如果是首层或 Fried 半径变了，就重算协方差拉普拉斯与谱模态
        if idxscr==1 || r0scrn(idxscr)~=r0scrn(idxscr-1)
            % 协方差拉普拉斯
            L_G = constructCovarianceLaplacian(pts, r0scrn(idxscr), alpha);
            % 取最小特征值对应的模式（长尺度模态）
            [PsiG, SigmaG] = eigs(L_G, M_turb, 'smallestabs', opts);
        end
        % KL 随机系数，复数形式
        sigma_n = sqrt(diag(SigmaG)/2);
        s_n = sigma_n .* (randn(M_turb,1) + 1i*randn(M_turb,1));
        phi = real(PsiG * s_n);

        % —— 幅度按结构函数缩放 —— 
        % 论文结构函数 D_phi(r)=6.88*(r/r0)^(5/3)
        % 对应相位谱强度按 (6.88)^(1/2)*(dz/r0)^(5/6)
        scale_factor = sqrt(6.88) * (dz / r0scrn(idxscr))^(5/6);
        phi = phi * scale_factor;

        %% —— 分裂步进 —— 
        % 前半步已经完成在循环外

        % 叠加湍流相位
        U_current = U_current .* exp(1i * phi);

        % 后半步衍射
        U_current = propagate_HalfStep(U_current, Phi, Lambda, dz/2, k);

        %% —— 自适应网格（每2层一次） —— 
        if mod(idxscr,2)==0
            % 重构网格
            [pts, TR] = adaptive_remeshing(pts, U_current, 30, 1, -30);
            % 重建传播图拉普拉斯与算子
            L_C = cotangent_Graph_Laplacian(pts, TR);
            [Phi, Lambda] = eigs(L_C, size(pts,1), 'smallestabs', opts);
            % 清除湍流模态缓存，以便下一层重算 L_G
            clear PsiG SigmaG L_G
        end

        % 保存每层后场（可视化或后处理）
        U_layers{idxreal, idxscr} = U_current;
    end

    % 最终场
    U_layers{idxreal, scr_count+1} = U_current;

    % 误差计算
    epsilon(idxreal) = compute_error(pts, TR, U_current);
end