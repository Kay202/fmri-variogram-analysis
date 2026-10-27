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
