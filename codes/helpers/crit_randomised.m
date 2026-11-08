% crit_randomised.m -- part of fmri-variogram-analysis
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

function [cstar, gam] = crit_randomised(mx, alpha)
%CRIT_RANDOMISED  Critical cluster size and randomisation probability for
%     a discrete max-cluster distribution at level alpha. The extent
%     distribution is discrete, so an exact test needs the randomisation
%     term as well as the threshold.
v = unique(mx);
cstar = v(end); gam = 1;
for i = numel(v):-1:1
    if mean(mx > v(i)) <= alpha && mean(mx >= v(i)) >= alpha
        cstar = v(i); break
    end
end
pgt = mean(mx >  cstar);
peq = mean(mx == cstar);
if peq > 0, gam = min(max((alpha - pgt)/peq, 0), 1); else, gam = 0; end
end

