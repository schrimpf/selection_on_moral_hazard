function r = mnrnd(n,p,m)
  persistent warned;
  if (~isOctave) % use Matlab builtin
    if (isempty(warned))
      warned = 0;
    end
    if (warned < 10) 
      warning(['You should remove mnrnd.m and use Matlab''s builtin ' ...
               'version']);
      warned = warned+1;
    end
  end
  
  if (nargin~=3) 
    error('3 inputs required');
  end
  if (n~=1)
    error('only n=1 implemented');
  end
  k = numel(p);
  if (size(p,2)~=k)
    p = p';
    if (size(p,2)~=k)
      error('p must be a vector');
    end
  end
  if (norm(sum(p)-1)>10*eps)
    error('sum of p is not 1');
  end
  r = zeros(m,k);
  u = rand(m,1)*ones(1,k);
  cp = ones(m,1)*cumsum(p);
  indj = sum(u > cp,2)+1;
  r(sub2ind(size(r),(1:m)',indj)) = 1;
end
