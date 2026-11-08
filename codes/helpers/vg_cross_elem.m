% vg_cross_elem.m -- part of fmri-variogram-analysis
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

function e = vg_cross_elem(S)
% On-axis cross structuring element of radius S.
%
% Eroding by this guarantees every pure-axis lag |k| <= S stays inside the mask,
% while retaining far more brain than a full (2S+1)^3 box -- about 63% against
% 34% at S=3 on these masks, which is why the decisive fixed-set test uses it.

e = false(2*S+1, 2*S+1, 2*S+1);
c = S + 1;
e(c, c, c) = true;
for k = 1:S
    e(c+k, c, c) = true;  e(c-k, c, c) = true;
    e(c, c+k, c) = true;  e(c, c-k, c) = true;
    e(c, c, c+k) = true;  e(c, c, c-k) = true;
end
end
