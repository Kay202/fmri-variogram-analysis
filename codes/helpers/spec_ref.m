% spec_ref.m -- part of fmri-variogram-analysis
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

function S = spec_ref(K, s, a, tau, sill)
% Closed-form spectral density of sum_m s_m Omega_3(a_m r) exp(-r^2/(2 tau_m^2)),
% evaluated at every mode. Non-negative by construction because (a-k)^2 < (a+k)^2
% for positive a and k. Only the shape matters: the scale is set from the sill,
% exactly as spec_gauss and spec_exp do.
Kp = max(K, eps);
P = zeros(size(K));
for m = 1:numel(s)
    P = P + s(m)*tau(m)*( exp(-((a(m)-Kp).^2)*tau(m)^2/2) ...
                        - exp(-((a(m)+Kp).^2)*tau(m)^2/2) ) ./ (Kp*a(m));
end
P(1) = 0;                 % the constant mode only shifts the field
P(P < 0) = 0;             % guards rounding only; the formula cannot be negative
S = sqrt(P * sill / sum(P(:)));
end

