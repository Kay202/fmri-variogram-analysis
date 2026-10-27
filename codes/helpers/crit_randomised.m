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

