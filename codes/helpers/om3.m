% om3.m -- part of fmri-variogram-analysis
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

function o = om3(x)
%OM3  Omega_3(x) = sin(x)/x, the isotropic basis function in three
%     dimensions, with the removable singularity at zero filled in.
o = ones(size(x)); nz = x ~= 0; o(nz) = sin(x(nz))./x(nz);
end

