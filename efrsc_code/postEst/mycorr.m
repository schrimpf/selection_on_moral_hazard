function c = mycorr(x,y);
  if (numel(x)==numel(y))
    c = corr(x(:),y(:));
  elseif size(x,2)==numel(y)
    Y = ones(size(x,1),1)*y(:)';
    c = corr(x(:),Y(:));
  end
end