%% 误差计算函数
function epsilon = compute_error(pts, TR, U)
% COMPUTE_ERROR 计算自适应网格下的仿真误差指标
%   epsilon = compute_error(pts, TR, U)
%   输入：
%       pts - N×2 顶点坐标
%       TR  - delaunayTriangulation 对象
%       U   - N×1 复振幅波前
%   输出：
%       epsilon - 标量，误差估计 ϵ = (1/N)∑_i h_i^2 |∇²U|_i

    % 1) 构造余切图拉普拉斯矩阵 L_C
    L_C = cotangent_Graph_Laplacian(pts, TR);
    % 2) 离散 Laplacian 应用到波前 U 上，近似 ∇²U
    LU = L_C * U;

    % 3) 计算每个顶点的平均边长 h_i
    E = TR.edges;  % M×2
    % 各边长度
    d = pts(E(:,1),:) - pts(E(:,2),:);
    len = sqrt(sum(d.^2, 2));
    % 累加到顶点
    N = size(pts,1);
    sum_len = accumarray([E(:,1);E(:,2)], [len;len], [N,1]);
    deg     = accumarray([E(:,1);E(:,2)], 1,     [N,1]);
    h_i     = sum_len ./ deg;

    % 4) 按公式计算误差 
    %    ε ≈ (1/N) ∑ h_i^2 * |(L_C U)_i|
    epsilon = mean( (h_i.^2) .* abs(LU) );
end
