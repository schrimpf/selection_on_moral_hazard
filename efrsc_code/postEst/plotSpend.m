function [q fy]=plotSpend(x,y,nq,fn,levels,makeFig,bw)
  if (nargin<5)
    levels = false;
  end
  if (nargin<6)
    makeFig = false;
  end
  if (nargin<7)
    bw = -1;
  end
  color = 'rbgkrbgkrbgk';
  lt = {'-','--','-.',':','-','--','-.',':','-','--','-.',':'};
  if (size(y,1)~= size(x,1))
    y = y';
  end
  assert(size(y,1)== size(x,1));
  for c=1:size(y,2)
    if (size(x,2)>1)
      [fy(:,c) q(:,c) bx(:,c)] = conditionalFoo(x(:,c),y(:,c),fn,nq,bw);
    else
      [fy(:,c) q(:,c) bx(:,c)] = conditionalFoo(x,y(:,c),fn,nq,bw);
    end
    if (nargout==0 || makeFig)
      if levels
        x = (bx(2:end,c)+bx(1:(end-1),c))/2;
        x(1) = bx(2,c)-1;
        line(x,fy(:,c),'Color',color(c),'LineStyle',lt{c});
      else
        line(q(:,c),fy(:,c),'Color',color(c),'LineStyle',lt{c});
      end
    end
  end
end

function [fy q bx]=conditionalFoo(x,y,foo,nq, bw)
% returns fy = foo(y|quantile(x)=q)
% for quantiles {k/nq}
  if (nargin<4)
    nq = 25;
  end
  if (nargin<3)
    foo = @(x) mean(x);
  end
  b = 1/nq:(1/nq):1;
  q = b-1/(2*nq);
  if (nargin<5 || bw<=0)
    b1x = quantile(x,b);
    b1x = b1x(:);
    b0x = [-Inf; b1x(1:end-1)];
  else
    b1x = quantile(x,min(q+bw/2,1));
    b1x = b1x(:);
    b0x = quantile(x,max(q-bw/2,0));
    b0x = b0x(:);
  end
  % brackets around each quantile
  fy = zeros(size(q));
  for j=1:numel(q)
    if (any(x>b0x(j) & x<=b1x(j)))
      fy(j) = foo(y(x>b0x(j) & x<=b1x(j)));
    else
      fy(j) = NaN;
    end
  end
  bx = [b0x(:); b1x(end)];
end
