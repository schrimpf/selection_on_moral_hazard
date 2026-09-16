function c = mycov(x,y);
  if (numel(x)==numel(y))
    C = cov(x,y);
    c = C(1,2);
  else
    c = 0;
    T = size(x,1);
    for t=1:T
      C = cov(x(t,:)',y);
      c = c+C(1,2);
    end
    c = c/T;
  end
end