%% 分谐波强化的FT相位屏生成
function [phz_lo, phz_hi] = ft_sh_phase_screen(r0,N,delta,L0,l0)
    D = N * delta;
    phz_hi = ft_phase_screen(r0,N,delta,L0,l0);
    % 空间坐标
    [x, y] = meshgrid((-N/2:N/2-1) * delta);
    
    phz_lo = zeros(size(phz_hi));
    for p = 1:3
        % 子谐波的频率步长
        del_f = 1/(3^p * D);
        fx=(-1:1)*del_f;
        [kx1,ky1]=meshgrid(fx);
        [th,f]=cart2pol(kx1,ky1);
        % 计算对应的 PSD 并生成 9 个复数系数
        fm = 5.92 / l0/(2 * pi);
        f0 = 1 / L0;
        PSD_phi = 0.023 * r0^(-5/3) .* exp(-(f/fm).^2)./(f.^2+f0.^2).^(11/6);
        PSD_phi(2,2) = 0;  % 避免零频发散
        cn = (randn(3)+1i*randn(3)).*sqrt(PSD_phi) * del_f;          %似乎有问题
        % cn = (randn(3) + 1i*randn(3)) .* sqrt(PSD_phi) * (del_f^2);     %用这个直接变平
        SH = zeros(N);
        for ii = 1:9
            SH = SH+cn(ii)*exp(1i*2*pi*(kx1(ii)*x+ky1(ii)*y));
        end
        phz_lo = phz_lo + SH;
    end
    phz_lo = real(phz_lo - mean(real(phz_lo(:))));
end
