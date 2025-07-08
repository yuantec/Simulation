function U0 = initialize_Gaussian_Beam(pts, w0, x0, y0)
%INITIALIZEGAUSSIANBEAM 在不规则三角网格顶点上生成高斯光束初始复振幅
%   U0 = initializeGaussianBeam(pts, w0, x0, y0)
%   输入:
%       pts - N×2 顶点坐标 [x, y]
%       w0  - 光束腰半径 (如 0.05 m)
%       x0  - 光束中心在 x 方向的偏移 (可设为 0)
%       y0  - 光束中心在 y 方向的偏移 (可设为 0)
%   输出:
%       U0  - N×1 复振幅向量，对应每个顶点处的初始光场

    % 提取顶点坐标
    x = pts(:,1);
    y = pts(:,2);
    % 计算相对于光束中心的坐标
    xp = x - x0;
    yp = y - y0;
    % 高斯振幅: exp(- (x^2 + y^2) / w0^2 )
    amp = exp(-(xp.^2 + yp.^2) / (w0^2));
    % 初始相位设为零，如需倾斜相位可在此叠加
    phase = zeros(size(amp));
    % 构造复振幅 U0 = A * exp(i * phi)
    U0 = amp .* exp(1i * phase);
end
