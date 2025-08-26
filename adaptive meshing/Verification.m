clc; clear; close all;
% 生成网格
N = 1024;
R = 1.0;
[pts, TR] = generate_Fibonacci_mesh(N, R);

% 构造 L
L = cotangent_Graph_Laplacian(pts, TR);

% 先做轻微数值对称化以消除浮点噪声（可选）
L = (L + L')/2;

% 对称性
symErr = norm(full(L - L'), 'fro');   % 确保传入 double
fprintf('symmetry error (fro): %.3e\n', symErr);

% 行和接近 0
rowSums = sum(L,2);                    % 可能为稀疏
rowSumMax = full(max(abs(rowSums)));   % 转成 full
fprintf('max abs row sum: %.3e\n', rowSumMax);

% 常向量是否在零空间
z = L * ones(size(pts,1),1);           % 乘出来通常是 full，但保险起见：
z = full(z);
fprintf('norm(L*1): %.3e\n', norm(z));

% 断言（使用 full/double）
assert(symErr < 1e-8, 'L 非对称（超过容差）');
assert(rowSumMax < 1e-8, '行和不为零（超过容差）');
assert(norm(z) < 1e-6, '常向量不是零特征向量（超过容差）');

% 额外：随机 Rayleigh 商（也把结果 full）
x = randn(N,1);
x = x - mean(x);
rq = (x'*(L*x)) / (x'*x);
rq = full(rq);
fprintf('random Rayleigh quotient: %.3e\n', rq);
if rq < -1e-12
    warning('发现显著负的 Rayleigh 商（可能存在负特征值或数值问题）');
end
