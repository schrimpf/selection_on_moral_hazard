%% Gaussian Kernel Density Estimator
% [bandwidth, density, xmesh, cdf] = kde(data, n, MIN, MAX)
function [bandwidth, density, xmesh, cdf] = kde(data, n, MIN, MAX)
  data = data(:);
  data = data(~isnan(data) & ~isinf(data));
  N = length(data);
  if N == 0
    bandwidth = 1; density = zeros(1, n); xmesh = linspace(0, 1, n); cdf = zeros(1, n);
    return;
  end
  if nargin < 2 || isempty(n)
    n = 2^14;
  end
  if nargin < 3 || isempty(MIN)
    MIN = min(data);
  end
  if nargin < 4 || isempty(MAX)
    MAX = max(data);
  end
  
  % Silverman's rule of thumb bandwidth
  s = std(data);
  iqr_val = iqr(data);
  if iqr_val > 0
    h = 0.9 * min(s, iqr_val / 1.34) * N^(-0.2);
  else
    h = 1.06 * s * N^(-0.2);
  end
  if h <= 0
    h = 1.0;
  end
  bandwidth = h;
  xmesh = linspace(MIN, MAX, n);
  
  % Gaussian kernel density evaluation
  % density(x) = 1/(N*h) * sum( normpdf((x - data)/h) )
  density = zeros(1, n);
  for i = 1:n
    u = (xmesh(i) - data) / h;
    density(i) = sum(exp(-0.5 * u.^2)) / (N * h * sqrt(2 * pi));
  end
  if nargout > 3
    cdf = cumsum(density) * (xmesh(2) - xmesh(1));
  end
end
