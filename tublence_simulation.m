%% 独立执行：点源湍流仿真完整脚本 
clc; clear; close all; 
 
%% —— 基本物理与系统参数 —— 
D2 = 0.5;            % 观测孔径直径 (m) 
wvl = 1e-6;          % 波长 (m) 
k = 2*pi/wvl;        % 波数 
Dz = 50e3;           % 总传播距离 (m) 
R = Dz;              % 波前曲率半径 
 
%% —— 点源与采样格距 —— 
d1 = 10e-3;       % 源面采样间距 (m) 
d2  = 10e-3;      % 观测面采样间距 (m) 
N = 512;          % 网格点数 N×N 
delta1 = d1; 
deltan = d2; 
 
%% —— 湍流相关量 —— 
Cn2   = 1e-16;         % 湍流结构常数 
L0 = 100;              % 外尺度 
l0 = 0.01;             % 内尺度 
r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);   % 平面波 Fried 参量 
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);       % 球面波 Fried 参量 
p = linspace(0,Dz,1e3); 
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1))); 
 
%% —— 优化求解各相位屏 r0scrn —— 
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
 
%% —— 构造源面场 pt —— 
DROI = 4*D2; 
D1 = wvl*Dz/DROI; 
[x1,y1] = meshgrid((-N/2:N/2-1)*delta1); 
[~,r1] = cart2pol(x1,y1); 
pt = exp(-1i*k/(2*R)*r1.^2)/D1^2 .* sinc(x1/D1) .* sinc(y1/D1) .* exp(-(r1/(4*D1)).^2); 
 
 
%% —— 计算 z 与混合采样 delta —— 
sg = exp(-(x1/(0.47 * N * d1)).^16) .* exp(-(y1/(0.47 * N * d1)).^16); 
t = repmat(sg,[1,1,n]); 
 
z = (1:n-1) * Dz/(n-1); 
zt = [0,z]; 
delta_z = zt(2:n)-zt(1:n-1); 
alpha = zt/zt(n); 
delta = (1-alpha) * delta1+alpha * deltan; 
 
phz = zeros(N,N,n); 
nreals = 40; 
[xn,yn,~] = ang_spec_multi_prop(pt,wvl,delta1,deltan,z,t); 
 
Uout = zeros(N); 
mask = circ(xn/D2,yn/D2,1); 
MCF2 = zeros(N); 
 
sg = repmat(sg,[1,1,n]); 
for idxreal = 1:nreals 
    for idxscr = 1:1:n 
        [plz_lo,phz_hi] = ft_sh_phase_screen(r0scrn(idxscr),N,delta(idxscr),L0,l0); 
        phz(:,:,idxscr) = plz_lo + phz_hi; 
    end 
    [xn,yn,Uout] = ang_spec_multi_prop(pt,wvl,delta1,deltan,z,sg.*exp(1i * phz)); 
    Uout = Uout .*exp(-1i * pi/(wvl * R) * (xn.^2+yn.^2)); 
    MCF2 = MCF2+corr2_ft(Uout,Uout,mask,deltan); 
end 
MCDOC2 = abs(MCF2)/(MCF2(N/2+1,N/2+1)); 
 
%% 画图 
% 1. 计算辐照度（强度）和相位 
Ivac = abs(Uout).^2;    % 辐照度分布 W/m^2 
Phivac = angle(Uout);     % 相位分布，单位 rad 
 
% 2. 构造归一化坐标轴（2x/D2, 2y/D2 范围约[-5,5]） 
%    注意：xn, yn 单位是米，除以 D2/2 再乘 2 ⇒ 2*xn./D2 
x_norm = xn(1, :); 
y_norm = yn(:, 1); 
 
% (a) 绘制辐照度 
figure(1); clf; 
imagesc(x_norm,y_norm, Ivac); 
axis image;                   % 保持 xy 比例 
set(gca,'YDir','normal');     % y 轴从下向上 
colormap(gray); 
cbar = colorbar; 
cbar.Label.String = '辐照度 (W/m^2)'; 
% 根据最大值微调 caxis：示意图中上限约 400，可按需调整 
clim([0, max(Ivac(:))]); 
xlabel('2x/D_2'); 
ylabel('2y/D_2'); 
title('(a) 辐照度分布'); 
 
% (b) 绘制相位 
figure(2); clf; 
imagesc(x_norm,y_norm, Phivac); 
axis image; 
set(gca,'YDir','normal'); 
colormap(gray); 
cbar = colorbar; 
cbar.Label.String = '相位 (rad)'; 
% 将相位限制在 ±π 
clim([-2, 2]); 
xlabel('2x/D_2'); 
ylabel('2y/D_2'); 
title('(b) 相位分布'); 

%% 结构函数验证
testLayer = round(scr_count/2); 
r0_layer = r0scrn(testLayer); 
delta_layer = delta(testLayer); 

xn2 = linspace(-N/2, N/2-1, N) * delta_layer;   % 创建坐标网格
yn2 = xn2;
[X2, Y2] = meshgrid(xn2);                       % 生成二维网格坐标

structure_fcn = @(phz, N) (sum((phz(:,2:N)-phz(:,1)).^2, 1)...
            +sum((phz(:,N-1:-1:1)-phz(:,N)).^2, 1)...
            +sum((phz(2:N,:)-phz(1,:)).^2, 2)'...
            +sum((phz(N-1:-1:1,:)-phz(N,:)).^2, 2)')/(N-1)/4;
NUM = 100;
SF_sh = zeros(1,N-1);
SF_hi = SF_sh;

for nn = 1:NUM
    [phz_lo, phz_hi] = ft_sh_phase_screen(r0_layer, N, delta_layer, L0, l0);
    phz = phz_lo + phz_hi;

    SF_sh = SF_sh + structure_fcn(phz, N);
    SF_hi = SF_hi + structure_fcn(phz_hi,N);
end
figure(3); clf; 
D_r = (0:N-2)*delta_layer;
plot(D_r/r0_layer,SF_sh/NUM)
hold on
plot(D_r/r0_layer,SF_hi/NUM)
plot(D_r/r0_layer, 6.88*(D_r/r0_layer).^(5/3))
xlim([0,12])
ylim([0,500])
legend('加次谐波','无次谐波','理论值')
xlabel('$|\Delta\textbf{r}|/{\it{r}}_0$', 'Interpreter', 'latex');
ylabel('$D_{\phi}(\Delta\textbf{r}) [\textbf{rad}^2]$', 'Interpreter', 'latex');

%% 相干因子验证
% 
% maxLag = floor(N/2);
% delta_layer = delta(end);     % 最后观测面采样间隔
% % 构建物理坐标矩阵
% coords = ((-N/2:N/2-1)*delta_layer);
% [X, Y] = meshgrid(coords, coords);
% R = sqrt(X.^2 + Y.^2);
% 
% % 径向平均
% emp_rho = 0:maxLag;
% emp_gamma = zeros(size(emp_rho));
% for d = emp_rho
%   ring = (R >= d*delta_layer) & (R < (d+1)*delta_layer) & mask;
%   emp_gamma(d+1) = mean(MCDOC2(ring),'all');
% end
% rho = emp_rho * delta_layer / r0sw;    % 归一化滞后
% 
% % 理论相干因子
% gamma_th = exp(-3.44 * (rho).^ (5/3));  % 或根据你的点源用合适系数
% 
% figure(4); clf;
% plot(rho, emp_gamma, 'bo-', 'DisplayName','仿真');
% hold on;
% plot(rho, gamma_th,'r-','LineWidth',1.5,'DisplayName','理论');
% hold off;
% xlabel('r/r_0');
% ylabel('相干因子 \gamma');
% legend('Location','northeast');
% grid on;
% xlim([0,1.5]);
% title('观测平面中相干因子：理论 vs 仿真');

% % 使用径向分箱方法计算相干因子
% ic = N/2+1; jc = N/2+1;  % 中心点坐标
% 
% % 计算每个点到中心点的距离（单位：米）
% dist_mat = sqrt((xn - xn(ic,jc)).^2 + (yn - yn(ic,jc)).^2);
% 
% % 设置径向分箱参数
% maxR = D2/2;          % 最大分析距离（孔径半径）
% dr = deltan;           % 径向分箱步长（与采样间距相同）
% rbins = 0:dr:maxR;     % 径向分箱边界
% rbins_center = (rbins(1:end-1) + rbins(2:end))/2;  % 分箱中心位置
% num_bins = length(rbins)-1;
% 
% % 初始化相干因子存储数组
% coherence_factor = zeros(1, num_bins);
% count = zeros(1, num_bins);
% 
% % 遍历所有有效点（掩膜内）
% for i = 1:N
%     for j = 1:N
%         if mask(i,j) == 1  % 只在有效区域内计算
%             r = dist_mat(i,j);
%             if r <= maxR
%                 % 确定所属径向分箱
%                 bin_idx = floor(r/dr) + 1;
%                 if bin_idx <= num_bins
%                     coherence_factor(bin_idx) = coherence_factor(bin_idx) + MCDOC2(i,j);
%                     count(bin_idx) = count(bin_idx) + 1;
%                 end
%             end
%         end
%     end
% end
% 
% % 计算平均相干因子
% coherence_factor = coherence_factor ./ max(count, 1);  % 避免除以零
% 
% % 理论相干因子模型 (球面波传播)
% gamma_theory_radial = exp(-3.44 * (rbins_center / r0sw).^(5/3));
% 
% % 绘制相干因子
% figure('Position', [100 100 800 600], 'Color', 'w', 'Name', '径向平均相干因子');
% plot(rbins_center/r0sw, coherence_factor, 'b--', 'LineWidth', 2, 'DisplayName', '仿真');
% hold on;
% plot(rbins_center/r0sw, gamma_theory_radial, 'r-', 'LineWidth', 2, 'DisplayName', '理论');
% hold off;
% 
% xlim([0 1.5]);
% ylim([0 1.0]);
% 
% % 添加标签和标题
% xlabel('\rho / r_0', 'FontSize', 14, 'Interpreter', 'tex');
% ylabel('相干因子', 'FontSize', 14);
% title('观察平面中的相干因子（径向平均）', 'FontSize', 16);
% 
% % 添加图例
% legend('Location', 'northeast', 'FontSize', 12);
% 
% 
% 
% % 添加网格
% grid on;
% box on;