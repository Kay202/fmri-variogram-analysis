function g = gen_gamma(Pw, nugget, V, k2u, gridsz)
%GEN_GAMMA  Exact variogram implied by a spectral power field, by inverse
%     FFT and averaging over the displacement vectors in each lag class.
C = real(ifftn(Pw)) * prod(gridsz);
g = zeros(numel(k2u),1);
for j = 1:numel(k2u)
    sel = V(V(:,4) == k2u(j), 1:3); acc = 0;
    for m = 1:size(sel,1)
        idx = mod(sel(m,:), gridsz) + 1;
        acc = acc + C(idx(1), idx(2), idx(3));
    end
    g(j) = nugget + (C(1,1,1) - acc/size(sel,1));
end
end

