function [ecost espend eoop] = exCost(data,nmh, negLambdaZeroUtil)
  
  if (nargin<2)
    nmh = false;
  end
  if (nargin<3)
    negLambdaZeroUtil = false;
  end

  [order xint wint] = gqzero(data.order);
  xint = xint'*sqrt(2);
  order = numel(xint);
  wint = wint/sqrt(pi);
  
  ecost = zeros(5,data.T,data.N);
  espend = ecost;
  eoop = ecost;
  td = data;
  for c=1:5
    td.choice(:) = c;
    for i=1:numel(xint);
      td.lambda = exp(td.muL' + ones(data.T,1)*td.sigL'*xint(i)) - ones(data.T,1)*td.lamlo';
      [totspend toop oopRate cost] = spending(td,nmh);
      ecost(c,:,:) = squeeze(ecost(c,:,:))  + cost * wint(i);
      espend(c,:,:) = squeeze(espend(c,:,:))  + totspend * wint(i);
      eoop(c,:,:) = squeeze(eoop(c,:,:)) + toop*wint(i);
    end
  end
  
end