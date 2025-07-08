N = 256;
L = 16;
delta = L/N;
F = 1/L;
x = (-N/2:N/2-1) * delta;
[x,y] = meshgrid(x);
w = 2;
A = rect(x/w) .* rect(y/w);
mask = ones(N);

C = str_fcn2_ft(A,mask,delta);
C_cont = 2 * w2 * (1-tri(x/w) .*tri(y/w));
