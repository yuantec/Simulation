%% 自适应网格重构函数
function [newPts, newTR] = adaptive_remeshing(pts, U, theta_max, theta_min, K_dB)
% ADAPTIVE_REMESHING 基于波前曲率对三角网格进行自适应重构
%   [newPts, newTR] = adaptive_remeshing(pts, U, theta_max, theta_min, K_dB)
%   输入：
%       pts       - N×2 网格顶点坐标 [x,y]\%       U         - N×1 波前复振幅，用于计算模值
%       theta_max - 曲率阈值，上限（度）
%       theta_min - 曲率阈值，下限（度）
%       K_dB      - 功率密度阈值（dB），用于顶点删除
%   输出：
%       newPts    - 重构后顶点坐标
%       newTR     - delaunayTriangulation 对象

    % 初始三角剖分
    TR = delaunayTriangulation(pts);
    T  = TR.ConnectivityList;
    N0 = size(pts,1);

    % 计算顶点法向量
    Vnorm = zeros(N0,2);
    for i = 1:N0
        triIdx = any(T==i,2);
        tris   = T(triIdx,:);
        verts  = pts(tris,:);
        % 面法向量2D: 取三角质心-顶点向量并旋转90°
        for j = 1:size(tris,1)
            v1 = verts(j,1:2) - pts(i,:);
            n1 = [ -v1(2), v1(1) ];
            Vnorm(i,:) = Vnorm(i,:) + n1;
        end
        % 单位化
        Vnorm(i,:) = Vnorm(i,:) / norm(Vnorm(i,:));
    end

    % 重构顶点集
    newPts = pts;
    for e = 1:size(TR.edges,1)
        a = TR.edges(e,1);
        b = TR.edges(e,2);
        % 计算两端法向量夹角
        cosTh = dot(Vnorm(a,:), Vnorm(b,:));
        theta = acosd(max(min(cosTh,1),-1));
        % 拆分或合并
        if theta > theta_max
            % 拆分：插入中点
            mid = mean(newPts([a,b],:),1);
            newPts(end+1,:) = mid;
        elseif theta < theta_min
            % 合并：删除邻接度低顶点
            % 延后删除，标记 a/b
            % 简化：只删除 b 点
            newPts(b,:) = [NaN NaN];
        end
    end
    % 删除标记顶点
    newPts = newPts(~any(isnan(newPts),2),:);

    % 根据功率密度删除低影响顶点
    Pd   = abs(U).^2;
    Pmean= mean(Pd);
    K    = 10^(K_dB/10);
    mask = Pd >= K*Pmean;
    newPts = newPts(mask,:);

    % 重新三角剖分
    newTR = delaunayTriangulation(newPts);
end
