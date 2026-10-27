function mx = max_clusters(S, nugget, gridsz, nsim, zthr, seed)
%MAX_CLUSTERS  Maximum suprathreshold cluster size over nsim Gaussian random
%     fields drawn from the given spectral amplitudes, at 26-connectivity.
rng(seed);
N = prod(gridsz);
mx = zeros(nsim, 1);
for i = 1:nsim
    f = real(ifftn(S .* (randn(gridsz) + 1i*randn(gridsz)))) * N;
    if nugget > 0, f = f + randn(gridsz) * sqrt(nugget); end
    z = (f - mean(f(:))) / std(f(:));
    cc = bwconncomp(z > zthr, 26);
    if cc.NumObjects == 0
        mx(i) = 0;
    else
        mx(i) = max(cellfun(@numel, cc.PixelIdxList));
    end
end
end
