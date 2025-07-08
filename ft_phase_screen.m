    %% 利用ft方法生成相位屏
    function phz = ft_phase_screen(r0,N,delta,L0,l0)
        del_f=1/(N*delta);                        % 频率间隔
        fx=(-N/2:N/2-1)*del_f;                    % 横向频率
        [kx ,ky]=meshgrid(fx);                    % 生成网格坐标矩阵
        [th ,f]=cart2pol(kx, ky);                 % 直极转化
        fm = 5.92/l0/(2 * pi);                    % 内尺度频率 [1/m]
        f0 = 1/L0;                                % 外尺度频率 [1/m]
        PSD_phi = 0.023 * r0^(-5/3).*exp(-(f/fm).^2)./(f.^2+f0.^2).^(11/6);
        PSD_phi(N/2+1,N/2+1) = 0;
        % 合成频谱系数
        cn = (randn(N) + 1i*randn(N)) .* sqrt(PSD_phi) * del_f;
        % 相位屏生成
        phz = real(ift2(cn,1));                  %可能是错的，所以改了一下
        % phz = real(ift2(cn, del_f));           %这个用了误差大   

