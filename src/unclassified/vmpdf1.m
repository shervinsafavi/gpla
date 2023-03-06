function y = vmpdf1(theta, mu, kappa)
% y = vmpdf1(theta, mu, sigma)

y = exp(kappa * cos(theta - mu)) / (2 * pi * besseli(0, kappa));
