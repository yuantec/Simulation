function y = tri(t)
    t = abs(x);
    y = zeros(size(t));
    idx = find(t<1.0);
    y(idx) = 1.0-t(idx);