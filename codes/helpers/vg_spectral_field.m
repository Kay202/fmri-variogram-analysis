function [S, info] = vg_spectral_field(p, nodes, gridsz, vox)
% Build a spectral amplitude mask for fast isotropic field synthesis.
%
% E1/E10 needs thousands of null realisations, and the random-wave synthesiser
% (vg_synth.m) costs ~20 s each -- too slow. This does it with one FFT.
%
% The Shapiro-Botha fit is already a spectral model:
%
%   C(r) = p0*1{r=0} + sum_k p_k * Omega_3(t_k r),   Omega_3(x) = sin(x)/x
%
% and Omega_3(t|h|) is the covariance of a plane wave of angular frequency t with
% uniformly random direction -- i.e. the spectral measure is a sum of non-negative
% delta shells at |k| = t_k. So the field can be drawn by putting power p_k on the
% Fourier modes lying on each shell, multiplying by complex Gaussian noise and
% inverse-transforming.
%
% This is exact for this model family, unlike circulant embedding of C(r), which
% fails here because sin(tr)/(tr) decays as 1/r and the wrapped covariance is not
% non-negative definite (§6.10).
%
% Non-negativity of p is what makes it legitimate: power cannot be negative, so a
% valid variogram is precisely what permits this construction.
%
% Returns S, the per-mode amplitude, to be used as
%   f = real(ifftn(S .* (randn(sz) + 1i*randn(sz)))) * info.scale
%
% with info.scale = prod(gridsz). MATLAB's ifftn divides by N, so for
% A_k = sqrt(P_k)(a+ib) the real part has variance sum(P)/N^2; multiplying by N
% recovers variance sum(P), which is the model's sill.

if nargin < 3 || isempty(gridsz), gridsz = [64 64 64]; end
if nargin < 4 || isempty(vox),    vox    = 3;          end

w = p(2:end); w = w(:);
nodes = nodes(:);
keep = w > 0;
w = w(keep); tk = nodes(keep);

% angular frequency of each Fourier mode, rad/mm
kk = cell(1,3);
for i = 1:3
    n = gridsz(i);
    f = [0:floor(n/2), -(ceil(n/2)-1):-1] / (n * vox);   % cycles/mm
    kk{i} = 2 * pi * f;
end
[K1, K2, K3] = ndgrid(kk{1}, kk{2}, kk{3});
Kmag = sqrt(K1.^2 + K2.^2 + K3.^2);

% assign each node's mass to the modes nearest its shell radius
P = zeros(gridsz);
edges = shell_edges(tk);
for j = 1:numel(tk)
    sel = Kmag >= edges(j) & Kmag < edges(j+1);
    n = nnz(sel);
    if n == 0
        [~, ix] = min(abs(Kmag(:) - tk(j)));   % nearest single mode
        sel = false(gridsz); sel(ix) = true; n = 1;
    end
    P(sel) = P(sel) + w(j) / n;
end

% The DC mode only shifts the field by a constant and contributes nothing to the
% variogram; dropping it avoids a spurious offset.
kmin = min(Kmag(Kmag > 0));
lost = P(1,1,1);
P(1,1,1) = 0;

S = sqrt(P);
info = struct('n_nodes', numel(tk), 'total_mass', sum(w), 'nugget', p(1), ...
              'modes_used', nnz(P > 0), 'grid', gridsz, 'vox', vox, ...
              'scale', prod(gridsz), ...
              'k_fundamental', kmin, ...
              'mass_below_fundamental', lost, ...
              'frac_mass_lost', lost / sum(w));
end

function e = shell_edges(t)
% midpoints between consecutive node radii, so every mode is claimed once
t = sort(t(:));
if isscalar(t)
    e = [0; inf]; return;
end
mid = (t(1:end-1) + t(2:end)) / 2;
e = [0; mid; inf];
end
