function dp = findcfPrice(eCost, dprem,subsidy, w)
% inputs: dprem = 1 x N , ecost = 2 x N, w = size(dprem)
  if (nargin<3)
    subsidy=0;
  end
  if (size(dprem,1)~=1)
    dprem = dprem';
  end
  if (nargin<4) 
    w = ones(size(dprem));
  end
  if (size(w,1)~=1);
    w = w';
  end
  valid = isfinite(dprem) & isfinite(eCost(1,:)) & isfinite(eCost(2,:)) & isfinite(w);
  dprem = dprem(valid);
  eCost = eCost(:, valid);
  w = w(valid);
  ec1 = @(dp) mean(eCost(1,dprem<=dp*(1-subsidy)).*w(dprem<=dp*(1-subsidy))) ...
        / mean(w(dprem<=dp*(1-subsidy)));
  ec5 = @(dp) mean(eCost(2,dprem<=dp*(1-subsidy)).*w(dprem<=dp*(1-subsidy))) ...
        / mean(w(dprem<=dp*(1-subsidy)));
  obj = @(dp) ec5(dp) - (ec1(dp) - dp);
  
  % find zero of obj. first construct initial bracket
  tau0 = 0.5;
  dp0 = quantilew(dprem',tau0,w)/(1-subsidy);
  f0 = obj(dp0);
  if isnan(f0)
   error('this is impossible');
  end
  dt = 0.005;
  tau1 = tau0 + dt;
  dp1 = quantilew(dprem',tau1,w)/(1-subsidy);
  iter = 0;
  while (~(obj(dp1)*f0<0))
    if (mod(iter,2)==1)
      dt = -dt;
    else 
      dt = abs(dt) + 0.005;
    end    
    iter = iter+1;
    tau1 = tau1+dt;
    if (tau1>=1 || tau1<=0) 
      fprintf('warning: tau1 out of range, failed to find bracket\n');
      if (f0>0) 
        fprintf('ecost(5) is always greater than ecost(1)\n');
        fprintf('setting dp such that all choose 1 (this dp is not unique)\n');
        dp = (min(dprem)-1)/(1-subsidy);
      else
        fprintf('ecost(1) is always greater than ecost(5)\n');
        fprintf('setting dp such that all choose 5 (this dp is not unique)\n');
        dp = (max(dprem)+1)/(1-subsidy);
      end
      return;
    end
    dp1 = quantilew(dprem',tau1,w)/(1-subsidy);
  end
  % find zero
  dp = fzero(obj, [dp0 dp1]);
end

%% quantile function with optional weights
function q=quantilew(x,tau,w) 
  if (nargin<3 || isempty(w) || isequal(w,ones(size(w))))
    q = quantile(x,tau);
    return;
  end
  
  [xs is] = sort(x);
  ws = w(is);
  cws = cumsum(ws)/sum(ws);
  
  q = tau;
  for i=1:numel(tau)
    hi = find(cws>tau(i),1);
    if isempty(hi)
      q(i) = xs(end);
    elseif hi==1
      q(i) = xs(hi);
    else % linearly interpolate
      q(i) = (xs(hi)*(cws(hi)-tau(i)) + xs(hi-1)*(tau(i)-cws(hi-1))) ...
             / (cws(hi)-cws(hi-1));
    end
  end
end