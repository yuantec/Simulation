%% 建立点源和接收平面结构关系和湍流相关量的MATLAB编码
D2 = 0.5;                   % 观测孔径直径
wvl = 1e-6;                 % 光的波长
k = 2 * pi/wvl;              % 波数
Dz = 50e3;                  % 传播距离

%% 用 sinc 模拟点光源
DROI=4 * D2;                  % 观测平面上​​感兴趣区域（4是为了覆盖主要光强分布）
D1 = wvl * Dz/DROI;           % 光场​​中心瓣的宽度
R = Dz;                       % 波前的​​曲率半径（光源在位置 z=0，则在传播距离 z=Dz 处，波前的曲率中心仍在光源位置）

%% 湍流
Cn2 = 1e-16;                  %湍流结构常数，表示湍流强度
r0sw = (0.423 * k^2 * Cn2 * 3/8 * Dz)^(-3/5);           % 平面波r0
r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3/5);                 % 球面波r0
p = linspace(0,Dz,1e3);
rytov = 0.563 * k^(7/6) * sum(Cn2 * (1-p/Dz).^(5/6) .* p.^(5/6) * (p(2)-p(1)));

%% 相位屏
nscr = 11;                  % 相位屏数量
A =zeros(2,nscr);
alpha = (0:nscr-1)/(nscr-1);
A(1,:) = alpha.^(5/3);
A(2,:) = (1-alpha).^(5/6) .* alpha.^(5/6);
b = [r0sw.^(-5/3);  rytov/1.33 * (k/Dz)^(5/6)];

x0 = (nscr/3 * r0sw * ones(nscr,1)).^(-5/3);

fun = @(X) sum((A * X(:)-b).^2);

x1 = zeros(nscr,1);
rmax = 0.1;
x2 = rmax/1.33 *(k/Dz)^(5/6)./A(2,:);
x2(A(2,:)== 0) = 50^(-5/3);
[X,fval,exitflag,Output] = fmincon(fun,x0,[],[],[],[],x1,x2);

r0scrn = X.^(-3/5);
r0scrn(isinf(r0scrn)) = 1e6;

bp = A * X(:);[bp(1)^(-3/5)  bp(2) * 1.33 * (Dz/k)^(5/6)];
[r0sw rytov];
