% mode_freqs.m -- part of fmri-variogram-analysis
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

function K = mode_freqs(gridsz, vox)
%MODE_FREQS  Angular frequency magnitude |k| at every Fourier mode of a
%     grid with the given size and voxel spacing.
kk = cell(1,3);
for i = 1:3
    n = gridsz(i);
    f = [0:floor(n/2), -(ceil(n/2)-1):-1] / (n * vox);
    kk{i} = 2 * pi * f;
end
[K1,K2,K3] = ndgrid(kk{1}, kk{2}, kk{3});
K = sqrt(K1.^2 + K2.^2 + K3.^2);
end
