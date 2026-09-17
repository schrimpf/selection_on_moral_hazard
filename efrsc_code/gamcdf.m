%% Cumulative distribution function of Gamma distribution
% a = shape, b = scale (mean = a*b)
function p = gamcdf(x, a, b)
  if nargin < 3
    b = 1;
  end
  p = gammainc(x ./ b, a);
end
