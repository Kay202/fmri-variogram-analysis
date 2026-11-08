% max_clusters.m -- part of fmri-variogram-analysis
% Copyright (C) 2026 Godfred Arhin and Emmanuel Kubuafor
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.

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
