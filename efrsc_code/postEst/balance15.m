%% find prices such that ratio portion choose plan 5 when only select 1 and 5 offered.
function sdata=balance15(sdata,ratio)
  sdata.avail(2:4,2,:) = 0;
  sdata.deduct(sdata.avail~=1) = -999;
  [dp p1 p5 val]= pdiff(0,sdata);
  dprem = quantile(squeeze(val(1,2,:)-val(5,2,:))./sdata.psi,ratio);
  [dp p1 p5 val]= pdiff(dprem,sdata);
  fprintf('p1=%g p5=%g\n',p1,p5);
  sdata.prem(1,2,:) = sdata.prem(1,2,:) + dprem;
end

%% returns P(1)-P(5) when set prem(1)+dprem1
function [p p1 p5 val]=pdiff(dprem1,sdata)
  [order xint wint] = gqzero(sdata.order);
  xint = xint'*sqrt(2);
  order = numel(xint);
  wint = wint/sqrt(pi);
  prem = sdata.prem;
  prem(1,2,:) = prem(1,2,:) + dprem1;
  [choice val]= findChoices_c(log(sdata.omega),1, sdata.lamlo, 1, sdata.muL', ...
                              sdata.sigL, sdata.psi,xint,wint, ...
                              sdata.deduct,sdata.maxoop, prem, sdata.avail, ...
                              sdata.choice, sdata.totalSpend, ...
                              sdata.maxIter);
  %assert(sum(choice(2,:)==2 | choice(2,:)==3 | choice(2,:)==4)==0);
  p1 = mean(choice(2,:)==1);
  p5 = mean(choice(2,:)==5);
  p = p1-p5;
end