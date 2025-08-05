clc; clear; close all;
%% —— 基本物理与系统参数 ——  
D2 = 0.5;            % 观测孔径直径 (m)
wvl = 1e-6;           % 波长 (m)
k = 2*pi/wvl;       % 波数
Dz = 50e3;           % 总传播距离 (m)

%% —— 初始网格与采样设置 —— 
N = 512;                          % 顶点数量
d1 = 10e-3;                        % 源面采样间距 (m)
radius = (N*d1)/sqrt(pi);              % 圆盘半径 ≈ (N d1)/√π
[pts, TR] = generate_Fibonacci_mesh(N, radius);

% 可视化网格
figure(1);
triplot(TR);
axis equal;
title('初始 Fibonacci Delaunay 网格');
xlabel('x (m)'); ylabel('y (m)');

%% —— 普通高斯光束参数 —— 
z0 = 0;        % 初始面放在 z=0 
w0 = 0.02;     % 束腰 2 cm

% 计算初始光场（高斯光束）
U0 = Gaussian_Beam(pts, w0, z0, k);
I0 = abs(U0).^2;
figure(2);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I0, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar;
title('初始高斯光束辐照度分布');

%% —— 构建传播算子（图拉普拉斯） —— 
% 生成初始图拉普拉斯矩阵 
L_C = cotangent_Graph_Laplacian(pts, TR);
% 检查对称性与半正定性 
assert(issymmetric(L_C),'L_C 非对称');
assert(eigs(L_C,1,'smallestabs')>=-1e-10,'L_C 存在负特征值');
assert(all(abs(sum(L_C,2))<1e-10),'L_C 行和不为零');

% 特征分解 
opts_eigs = struct('issym',true,'isreal',false);
[Phi, Lambda] = eigs(L_C, N, 'smallestabs', opts_eigs);

%% —— 湍流相关物理量 & 分段 Fried 半径 —— 
alpha = 5/3;       % Kolmogorov 指数
Cn2 = 1e-16;     % 湍流结构常数

r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);   % 平面波 Fried 参量  
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);       % 球面波 Fried 参量 
p = linspace(0,Dz,1e3);
rytov= 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1)));

% 分段优化求解每层 r0scrn
scr_count = 11;               % 相位屏数量  
n = scr_count;
A = zeros(2,n);
alpha_vec = (0:n-1)/(n-1);
A(1,:) = alpha_vec.^(5/3);
A(2,:) = (1-alpha_vec).^(5/6).*alpha_vec.^(5/6);
b  = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
x0 = (n/3*r0sw*ones(n,1)).^(-5/3);
x1 = zeros(n,1);
rmax= 0.1;
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:);
x2(A(2,:)==0) = 50^(-5/3);

opts_fmin = optimoptions('fmincon','Display','none');
[X,~,~,~] = fmincon(@(X) sum((A*X - b).^2), x0,[],[],[],[], x1, x2,[], opts_fmin);
r0scrn = X.^(-3/5);
r0scrn(isinf(r0scrn)) = 1e6;

%% —— 多层分裂步进传播 —— 
M_turb = 100;                   % KL 模式数
dz  = Dz/scr_count;             % 每段距离
nreals  = 40;                   % 蒙特卡洛次数

% 初始半步衍射
U_half = propagate_HalfStep(U0, Phi, Lambda, dz/2, k);

% 预分配
U_layers = cell(nreals, scr_count+1);
epsilon  = zeros(nreals,1);

for real = 1:nreals
    U_current = U_half;
    U_layers{real,1} = U0; 
    
    for layer = 1:scr_count
        %—— 生成 KL 相位屏 —— 
        % 如果是首层或 Fried 半径变了，就重算协方差拉普拉斯与谱模态 
        if layer==1 || r0scrn(layer)~=r0scrn(layer-1)
            % 协方差拉普拉斯 
            L_G = construct_Covariance_Laplacian(pts, r0scrn(layer), alpha);
            % 取最小特征值对应的模式（长尺度模态）
            [PsiG, SigmaG] = eigs(L_G, M_turb, 'smallestabs', opts_eigs);
        end
        % KL 随机系数，复数形式 
        sigma_n = sqrt(diag(SigmaG)/2);
        s_n = sigma_n .* (randn(M_turb,1)+1i*randn(M_turb,1));
        phi = real(PsiG * s_n);

        % —— 幅度按结构函数缩放 ——  
        % 论文结构函数 D_phi(r)=6.88*(r/r0)^(5/3) 
        % 对应相位谱强度按 (6.88)^(1/2)*(dz/r0)^(5/6)
        scale_f = sqrt(6.88)*(dz/r0scrn(layer))^(5/6);
        phi = phi * scale_f;
        
        %—— 后半步衍射 & 相位叠加 —— 
        % 叠加湍流相位
        U_current = U_current .* exp(1i * phi);
        % 后半步衍射
        U_current = propagate_HalfStep(U_current, Phi, Lambda, dz/2, k);
        
        %—— 自适应网格（每2层一次） —— 
        if mod(layer,2)==0
            % 重构网格 
            [pts, TR] = adaptive_remeshing(pts, U_current, 30, -30, 1);
            % 重建传播图拉普拉斯与算子 
            L_C = cotangent_Graph_Laplacian(pts, TR);
            [Phi, Lambda] = eigs(L_C, size(pts,1), 'smallestabs', opts_eigs);
            % 清除湍流模态缓存，以便下一层重算 L_G 
            clear PsiG SigmaG L_G
        end

        % 保存每层后场（可视化或后处理）
        U_layers{real,layer} = U_current;
    end
    % 最终场 
    U_layers{real, scr_count+1} = U_current;
    % 误差计算
    epsilon(real) = compute_error(pts, TR, U_current);
end

%% 可视化示例：第1次试验的初始、中间和最终强度
for layer = [1, 6, scr_count+1]
    I = abs(U_layers{1,layer}).^2;
    figure;
    trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I, 'EdgeColor','none');
    view(2); shading interp; axis equal; colorbar;
    if layer==1
        title('试验1：初始场');
    elseif layer==scr_count+1
        title('试验1：最终场');
    else
        title(sprintf('试验1：第 %d 层后', layer-1));
    end
end
