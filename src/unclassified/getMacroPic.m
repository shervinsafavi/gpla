function [mu, kappa] = getMacroPic(svec);

% r = nanmean(svec, 1);
% aveLen = nanmean(abs(svec), 1)
% r = nanmean(svec / aveLen, 1);
r = nanmean(exp(1i * angle(svec)), 1);
cvar = 1 - abs(r);
mu = angle(r);
n = sum(~isnan(svec));
kappa = MLest_kappa_forVonMisesBasedPLV1(abs(r), n);

