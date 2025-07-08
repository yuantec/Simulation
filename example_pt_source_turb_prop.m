%% 通过采样分析给定网格的情况下，点源湍流仿真执行的MATLAB代码
l0 = 0;
L0 = inf;

zt = [0,z];
delta_z = zt(2:n)-zt(1:n-1);
alpha = zt/zt(n);
delta = (1-alpha) * delta1+alpha * deltan;

phz = zeros(N,N,n);
nreals = 20;

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