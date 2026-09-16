%IWISHRND Generate inverse Wishart random matrix
function [a,di] = iwishrnd(sigma,df,di)
  persistent warned;
  if (~isOctave) % use Matlab builtin
    if (isempty(warned))
      warned = 0;
    end
    if (warned < 10) 
      warning(['You should remove iwishrnd.m and use Matlab''s builtin ' ...
               'version']);
      warned = warned+1;
    end
  end
  
  if (nargin<3)
    if (isOctave)
      [d p] = chol_c(sigma);
    else 
      [d p] = chol(sigma);
    end
    if (p~=0)
      error('sigma not full rank');
    end
    di = d' \ eye(size(d,1));
  end
  x = randn(df,size(di,1))*di ;
  iA = x'*x;
  a = inv(iA);
end
