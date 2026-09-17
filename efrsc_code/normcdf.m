%% Normal cumulative distribution function
% p = normcdf(x, mu, sigma)
function p = normcdf(x, m, s)
  if nargin < 2, m = 0; end
  if nargin < 3, s = 1; end
  z = (x - m) ./ (s * sqrt(2));
  p = 0.5 * erfc(-z);
end
