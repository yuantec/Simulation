function L = cotangent_Graph_Laplacian(pts, TR)
% COTANGENT_GRAPH_LAPLACIAN 构造余切图拉普拉斯矩阵 L_C
%   L = cotangent_Graph_Laplacian(pts, TR)
%   输入：
%       pts - N×2 顶点坐标 [x, y]
%       TR  - delaunayTriangulation 对象
%   输出：
%       L   - N×N 稀疏余切拉普拉斯矩阵

    % 顶点数量 N 和三角形列表 T (M×3)
    N = size(pts,1);
    T = TR.ConnectivityList;
    M = size(T,1);  % 三角形数

    % 预估非零条目数: 每个三角形3条边，每条边2个非对角 + 2个对角
    maxEntries = M * 3 * 4;
    I = zeros(maxEntries,1);
    J = zeros(maxEntries,1);
    V = zeros(maxEntries,1);
    idx = 0;  % 当前填充位置

    % 遍历每个三角形
    for t = 1:M
        tri = T(t,:);
        % tri = [i, j, k]
        for e = 1:3
            % 边 (a,b) 及其对顶点 c
            a = tri(e);
            b = tri(mod(e,3)+1);
            c = tri(mod(e+1,3)+1);

            % 计算余切值 cot_alpha on triangle t
            u = pts(a,:) - pts(c,:);
            v = pts(b,:) - pts(c,:);
            cross_uv = cross([u,0],[v,0]);
            if norm(cross_uv) < eps
                continue;  % 跳过退化三角形导致的除零
            end
            cot_alpha = dot(u,v) / norm(cross_uv);

            % 查找共用边 (a,b) 的另一三角形 t2，得到第二个对顶点 d
            adjTris = edgeAttachments(TR, a, b);
            trisList = adjTris{1};
            % trisList 包含当前 t, 以及其它可能的 t2
            otherT = setdiff(trisList, t);
            if isempty(otherT)
                cot_beta = 0;  % 边界边，仅一个三角形时置零
            else
                % 取第一个相邻三角形
                t2 = otherT(1);
                tri2 = T(t2,:);
                d = tri2(~ismember(tri2, [a,b]));
                % 计算 cot_beta
                u2 = pts(a,:) - pts(d,:);
                v2 = pts(b,:) - pts(d,:);
                cross_u2v2 = cross([u2,0],[v2,0]);
                if norm(cross_u2v2) < eps
                    cot_beta = 0;
                else
                    cot_beta = dot(u2,v2) / norm(cross_u2v2);
                end
            end

            % 边权重 w = 0.5*(cot_alpha + cot_beta)
            w = 0.5 * (cot_alpha + cot_beta);

            % 将权重 w 累加到稀疏矩阵条目：
            % 非对角元 L(a,b)=L(b,a)=-w
            idx = idx + 1; I(idx)=a; J(idx)=b; V(idx)=-w;
            idx = idx + 1; I(idx)=b; J(idx)=a; V(idx)=-w;
            % 对角元   L(a,a)=L(a,a)+w, L(b,b)=L(b,b)+w
            idx = idx + 1; I(idx)=a; J(idx)=a; V(idx)= w;
            idx = idx + 1; I(idx)=b; J(idx)=b; V(idx)= w;
        end
    end

    % 裁剪多余预分配
    I = I(1:idx);
    J = J(1:idx);
    V = V(1:idx);

    % 构造稀疏拉普拉斯矩阵
    L = sparse(I, J, V, N, N);

    % 最后检查对称性和谱半正定性
    assert(issymmetric(L), '拉普拉斯矩阵非对称');
    minEig = eigs(L,1,'smallestabs');
    assert(minEig >= -1e-10, '拉普拉斯矩阵非半正定');
end
