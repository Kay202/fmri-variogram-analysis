function g = vg_axis_gamma(res2, sz, sel, kmax, tidx)
% Pure-axis semivariance on a fixed base voxel set, symmetrised over +k and -k.
%
%   res2  (T x nvox) single, time down the rows so the per-lag column gather
%         is contiguous under MATLAB's column-major layout
%   sz    [X Y Z] volume dimensions
%   sel   logical volume selecting the base voxels (must be eroded enough that
%         every x +/- k*e_a stays inside the mask for k <= kmax)
%   tidx  optional time indices, for the finite-sample analysis (E5); defaults
%         to all volumes
%   g     3 x kmax, rows are x/y/z
%
% Averaging +k and -k removes the first-order drift that arises because the
% partner set translates toward one side of the brain as the lag grows.

if nargin >= 5 && ~isempty(tidx)
    res2 = res2(tidx, :);
end
ci = find(sel);
[i, j, k] = ind2sub(sz, ci);
A = res2(:, ci);
g = nan(3, kmax);
for a = 1:3
    for m = 1:kmax
        acc = 0;
        for sgn = [1 -1]
            off = [0 0 0]; off(a) = sgn * m;
            si = sub2ind(sz, i + off(1), j + off(2), k + off(3));
            d = A - res2(:, si);
            acc = acc + 0.5 * mean(double(d(:)).^2);
        end
        g(a, m) = acc / 2;
    end
end
end
