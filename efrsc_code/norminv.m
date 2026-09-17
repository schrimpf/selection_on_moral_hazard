%% Inverse normal cumulative distribution function (probit)
% x = norminv(p, mu, sigma)
function x = norminv(p, m, s)
  if nargin < 2, m = 0; end
  if nargin < 3, s = 1; end
  x = m + s * sqrt(2) * erfinv(2 * p - 1);
end
