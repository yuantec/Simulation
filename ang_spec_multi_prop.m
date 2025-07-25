%% 菲涅尔衍射积分
function [xn,yn,Uout] = ang_spec_multi_prop(Uin,wvl,delta1,deltan,z,t)
    N = size(Uin,1);                   
    [nx,ny] = meshgrid(-N/2:N/2-1);    
    k = 2*pi/wvl;                      

    %% 超高斯吸收边界
    nsq = nx.^2 + ny.^2;
    w = 0.47 * N;
    sg = exp(-nsq.^8 / w^16);
    clear ('nsq','w');

    %% 构造分段距离
    z = [0, z];                        
    n = length(z);
    Delta_z = z(2:n) - z(1:n-1);

    %% 网格间距随距离线性变化
    alpha = z / z(n);                        % 1×n 向量
    delta = (1 - alpha) * delta1 + alpha * deltan;  
    m = delta(2:n)./delta(1:n-1);
    x1 = nx * delta(1);
    y1 = ny * delta(1);
    rlsq = x1 .^2 + y1.^2;
    Q1 = exp(1i*k/2 * (1-m(1))/Delta_z(1) * rlsq);
    
    %% ... 下面接着写每段传播的 Angular Spectrum 算法
    Uin = Uin .* Q1 .* t(:,:,1);
    for idx = 1:n-1
        % 计算 m 段传播用的频域坐标
        deltaf = 1/(N * delta(idx));
        fX = nx  * deltaf;
        fY = ny  * deltaf;
        fsq = fX.^2 + fY.^2;
        Z = Delta_z(idx);
        Q2 = exp(-1i * pi^2 * 2 * Z/m(idx)/k * fsq);
        Uin = sg .* t(:,:,idx+1) .* ift2(Q2.*ft2(Uin/m(idx),delta(idx)),deltaf);
    end

    % 计算输出坐标网格
    xn = nx * delta(n);
    yn = ny * delta(n);
    rnsq = xn.^2 + yn.^2;
    Q3 = exp(1i * k/2 * ((m(n-1)-1) / (m(n-1)*Z)) * rnsq);
    Uout = Q3 .* Uin;
end
