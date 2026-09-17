function plotP5psiOm(data,prefix)
  psi = quantile(data.psi,[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
  %psi = 10.^(-1:-1:-7);
  %lom = log(data.omega);
  %lo2 = (lom-mean(lom))*10 + mean(lom);
  %data.omega = exp(lo2);
  %data.muL = data.muL/2;
  %data.lamlo = data.lamlo + 1000;
  ind = data.choice(2,:)>0;
  for p=1:numel(psi)
    [c ev]= choicePsi(data,psi(p));
    choice(:,p) = (ev(5,2,ind) - ev(1,2,ind))/psi(p);
  end
  figure;
  plotSpend(data.omega(ind),choice,20,@(x) mean(x));
  xlabel('Quantiles of \omega');
  ylabel('P(c=5)');
  legend(sprintf('\\psi=%.2g',psi(1)), ...
         sprintf('\\psi=%.2g',psi(2)), ...
         sprintf('\\psi=%.2g',psi(3)), ...
         sprintf('\\psi=%.2g',psi(4)), ...
         sprintf('\\psi=%.2g',psi(5)), ...
         sprintf('\\psi=%.2g',psi(6)), ...
         sprintf('\\psi=%.2g',psi(7)));
  print('-depsc2',sprintf('figures/%sP5psiOm',prefix));
end

function [c ev]= choicePsi(data,psi)
  [order xint wint] = gqzero(data.order);
  xint = xint'*sqrt(2);
  order = numel(xint);
  wint = wint/sqrt(pi);

  dat = data;
  dat.psi = psi*ones(size(data.psi));
  [c ev es] = findChoices_c(log(dat.omega),1, dat.lamlo, 1, ...
                            dat.muL', dat.sigL, dat.psi, ...
                            xint,wint, dat.deduct, ...
                            dat.maxoop, dat.prem, dat.avail, ...
                            dat.choice, dat.totalSpend, ...
                            dat.maxIter);
end