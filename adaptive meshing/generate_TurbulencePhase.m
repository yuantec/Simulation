%% 生成湍流相位屏函数
function S = generate_TurbulencePhase(pts, rc, alpha, M)
% GENERATETURBULENCEPHASE 基于 Karhunen–Loève 展开生成湍流相位屏
%   S = generateTurbulencePhase(pts, rc, alpha, M)
%   输入:
%       pts   - N×2 网格顶点坐标 [x, y]
%       rc    - 波前相干半径
%       alpha - 相位相关函数指数 (通常 5/3)
%       M     - 谱模式数量 (M ≤ N)
%   输出:
%       S     - N×1 相位屏向量

    N = size(pts,1);
    % 1) 计算顶点两两距离
    Dmat = squareform(pdist(pts));

    % 2) 构造协方差权重 WΓ = Γs(r1,r2) = (|Δr|/rc)^alpha
    % 根据Kolmogorov模型用指数衰减近似协方差
    W = exp(-(Dmat / rc) .^ alpha);
    W(1:N+1:end) = 0;  % 对角置零

    % 3) 构造协方差拉普拉斯 L_G = Dg - W
    Dg = diag(sum(W,2));
    L_G = Dg - W;

    % 4) 特征分解 L_G = Ψ Σ Ψ^T
    opts.issym = true;
    [Psi, Sigma] = eigs(L_G, M, 'smallestabs', opts);

    % 5) Karhunen–Loève 随机系数 s_n ~ N(0, σ_n^2)
    sigma_n = sqrt(diag(Sigma));
    s_n = sigma_n .* randn(M,1);

    % 6) 构造相位屏 S = Ψ * s_n
    S = Psi * s_n;
end
