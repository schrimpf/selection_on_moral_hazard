%% Random draws from Gamma distribution
% Compatible with MATLAB gamrnd(a, b, varargin) where:
% a = shape, b = scale (mean = a*b)
function r = gamrnd(a, b, varargin)
  if nargin < 2
    error('gamrnd requires at least shape and scale parameters');
  end
  if nargin == 2
    sz = size(a);
    if isscalar(a) && ~isscalar(b)
      sz = size(b);
    end
    r = b .* randg(a, sz);
  else
    if length(varargin) == 1 && isvector(varargin{1}) && length(varargin{1}) > 1
      sz = varargin{1};
    else
      sz = [varargin{:}];
    end
    r = b .* randg(a, sz);
  end
end
