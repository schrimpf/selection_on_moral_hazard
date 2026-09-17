%% Normal probability density function
% y = normpdf(x, mu, sigma)
function y = normpdf(x, m, s)
  if nargin < 2, m = 0; end
  if nargin < 3, s = 1; end
  z = (x - m) ./ s;
  y = exp(-0.5 * z.^2) ./ (s * sqrt(2 * pi));
end
