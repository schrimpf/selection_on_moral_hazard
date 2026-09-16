%% returns spending -- data should be simulated data and include
%  data.lambda as a field
function [totspend toop oopRate cost spend util] = spending(data,nmh,negLambdaZeroUtil)

  if (nargin<2)
    nmh = false;
  end

  if (nargin<3)
    negLambdaZeroUtil = false;
  end

  T = data.T;
  N = data.N;
  ded = zeros(size(data.choice))*NaN;
  maxoop = ded;
  copay = ded;
  prem = ded;
  for c=1:5
    ded(data.choice==c) = data.deduct(c,data.choice==c);
    maxoop(data.choice==c) = data.maxoop(c,data.choice==c);
    copay(data.choice==c) = data.copay(c,data.choice==c);
    prem(data.choice==c) = data.prem(c,data.choice==c);
  end
  if (nmh)
    totspend = data.lambda;
    totspend(totspend<0) = 0;
    toop = totspend;
    oopRate = ones(size(toop));
    oopRate(totspend==0) = inf;
    oopRate(toop>ded) = 0.1;
    tmp = ded + 0.1*(totspend - ded);
    toop(toop>ded) = tmp(toop>ded);
    oopRate(toop>maxoop) = 0;
    toop(toop>maxoop) = maxoop(toop>maxoop);
    cost = totspend - toop - prem;
    return;
  end
  
  spend(1,:,:) = zeros(T,N);
  spend(2,:,:) = data.lambda;
  spend(3,:,:) = data.lambda+0.9*(data.omega*ones(1,T))';
  spend(4,:,:) = data.lambda+(data.omega*ones(1,T))';
  oop(1:2,:,:) = spend(1:2,:,:);
  oop(3,:,:) = ded + (squeeze(spend(3,:,:)) - ded)*0.1;
  oop(4,:,:) = maxoop;
  for j=1:4;
    L(j,:,:) = data.lambda;
    O(j,:,:) = (data.omega*ones(1,T))';
  end
  util = (spend-L)-(spend-L).^2./(2*O) - oop;
  util(2,squeeze(spend(2,:,:))<0 | squeeze(spend(2,:,:))>ded)=-Inf;
  util(3,squeeze(spend(3,:,:))<ded | squeeze(oop(3,:,:))>maxoop) = -Inf;
  util(4,squeeze(spend(4,:,:))<0 | (ded + (squeeze(spend(4,:,:)) - ...
                                           ded)*0.1<maxoop)) = -Inf;
  if (negLambdaZeroUtil) 
    util(1,data.lambda<0) = 0;
  end
  for j=1:4
    util(j,isnan(ded) | ded<0 | isnan(maxoop) | maxoop<0) = NaN;
  end
  [maxu seg]=max(util,[],1);
  totspend = zeros(size(data.lambda))*NaN;
  toop = zeros(size(data.lambda))*NaN;
  seg = squeeze(seg);
  seg(any(isnan(util),1))=NaN;
  for j=1:4
    totspend(seg==j) = spend(j,seg==j);
    toop(seg==j) = oop(j,seg==j);    
  end
  oopRate = inf*ones(size(data.lambda));
  oopRate(seg==2) = 1;
  oopRate(seg==3) = 0.1;
  oopRate(seg==4) = 0;
  cost = totspend - toop - prem;
end
