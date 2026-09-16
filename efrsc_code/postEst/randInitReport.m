clear;
close all;
% use parameters from estimates
config;
load(resultFile);
% compute posterior means
nSaved = numel(chain.theta);
post = floor(nSaved/2):nSaved;
for k=1:numel(chain.beta);
  beta0{k} = mean(chain.beta{k}(post,:),1)';
end
sig0 = squeeze(mean(chain.sig(post,:,:),1));
gamma0 = squeeze(mean(chain.gamma(post,:,:),1));
theta0 = squeeze(mean(chain.theta(post)));
shape0 = squeeze(mean(chain.shape(post)));
rho0 = squeeze(mean(chain.rho(post)));
%clear chain data;

names = {'\omega','\psi','\mu_{\lambda 1}','\mu_{\lambda 2}'};
xnames{1} = {'Constant','family','+ spouse','+ child', ...
             'switch 04','switch 05','switch 06','age', ...
             'female','tenure','income','rs\_fam','2004'};
xnames{2} = xnames{1};
xnames{3} = xnames{1};
%color = 'bgcmykbgcmyk';
lt =    '::::::::::::';
% 5 shades of orange
color{1}=[139,125,107]/255;
color{2}=[227,168,105]/255;
color{3}=[255,127,0]/255;
color{4}=[200,5,0]/255;
%color{5}=[255,128,0]/255;
% 5 shades of blue
color{5}=[106,90,205]/255;
color{6}=[61,89,171]/255;
color{7}=[72,118,255]/255;
color{8}=[162,181,205]/255;
%color{10}=[30,144,255]/255;
color{1} = 'b';

suffix = prefix;
%load('simTruncOm02'); %tmp_real0304only');
%data = simdata;
incTrue = false;

all{1} = chain;
%{
for j=1:8
  load(sprintf('simWithX%02d',j));
  all{j} = chain;
end

for j=1:5
  load(sprintf('oldmat/r1i%02d',j));
  all{j} = chain;
  load(sprintf('oldmat/r2i%02d',j));
  all{j+5} = chain;  
end
%}
clear chain;

t = (1:size(all{1}.beta{1},1))';
T = size(all{1}.beta{1},1);

first = 10;
t = first:T;
post = (floor(numel(t)/2):numel(t));
tau = [0.05 0.95];
for k=1:numel(all)
  fprintf('%2d: mean(sigOm)=%g\n',k,mean(all{k}.sig(floor(T/2):T,1,1)));
end
t = first:T;

sub1 = 3; sub2=2;
nsub = sub1*sub2;
sp = 1;
ns = 1;
for j=1:numel(all{1}.beta)
  bhi=size(all{1}.beta{j},2);
  for b=1:bhi
    fprintf('ihi=%d i=%d\n',bhi,b);
    for k=1:numel(all)    
      mu = all{k}.beta{j}(first:T,:);
      fprintf('j=%d i=%d k=%d\n',j,b,k);
      if (k==1) 
        if (sp==1)
          figure;
        end
        subplot(sub1,sub2,sp),plot(t,mu(:,b),'Color',color{k},'LineStyle', ...
                                   lt(k));
        sp = sp+1;
      else 
        line(t,mu(:,b),'LineStyle',lt(k),'Color',color{k});
      end
      avg = mean(mu(post,b));
      q = quantile(mu(post,b),tau);
      for ti=1:numel(tau) 
        line([first,T],[q(ti),q(ti)],'LineStyle','--','Color','k');
        text(first+(T-first)*7/10,q(ti),sprintf('Posterior %3.2f Quantile',tau(ti)));
      end
      line([first,T],[avg,avg],'LineStyle','-','Color','k');
      text(first + (T-first)*3/10,avg,'Posterior Mean');
    end % k
    if (incTrue)
      line([first,T],[beta0{j}(b),beta0{j}(b)],'LineStyle','-','Color','r')
      text(first+(T-first)/20,beta0{j}(b),'True Value');
    end
    title(sprintf('\\beta_%s(%s)',names{j},xnames{j}{b}));
    xlim([0 T]);
    %print('-depsc2',sprintf('figures/mu%d%s.eps',i,suffix));
    if (sp>nsub)
      print('-depsc2',sprintf('figures/%s%d',suffix,ns));
      ns = ns+1;
      sp = 1;
    end
  end % i
end % j

for i=1:size(all{1}.sig,2);
  for j=i:size(all{1}.sig,3);
    for k=1:numel(all)      
      sig = all{k}.sig(first:T,:,:);
      avg = mean(sig(post,i,j));
      q = quantile(sig(post,i,j),tau);
      if (k==1)
        if (sp==1) 
          figure;
        end
        subplot(sub1,sub2,sp),plot(t,sig(:,i,j),'Color',color{k},'LineStyle', ...
                                   lt(k));
        sp = sp+1;
      else 
        line(t,sig(:,i,j),'LineStyle',lt(k),'Color',color{k});
      end
      for ti=1:numel(tau) 
        line([first,T],[q(ti),q(ti)],'LineStyle','--','Color','k');
        text(first+(T-first)*7/10,q(ti),sprintf('Posterior %3.2f Quantile',tau(ti)));
      end
      line([first,T],[avg,avg],'LineStyle','-','Color','k');
      text(first+(T-first)*3/10,avg,'Posterior Mean');    
    end
    if (incTrue)
      line([first,T],[sig0(i,j),sig0(i,j)],'LineStyle','-','Color','r')
      text(first+(T-first)/20,sig0(i,j),'True Value');
    end
    title(sprintf('\\sigma_{%s,%s}',names{i},names{j}));
    xlim([0 T]);
    if (sp>nsub)
      print('-depsc2',sprintf('figures/%s%d',suffix,ns));
      ns = ns+1;
      sp = 1;
    end
    %print('-depsc2',sprintf('figures/sig%d%d%s.eps',i,j,suffix));
  end
end

for k=1:numel(all)
  rho = all{k}.rho(first:T);
  avg = mean(rho);
  q = quantile(rho(post),tau);
  if (k==1)
    if (sp==1) 
      figure;
    end
    subplot(sub1,sub2,sp),plot(t,rho,'Color',color{k},'LineStyle', ...
                               lt(k));
    sp=sp+1;
  else 
    line(t,rho,'LineStyle',lt(k),'Color',color{k});
  end
  for ti=1:numel(tau) 
    line([first,T],[q(ti),q(ti)],'LineStyle','--','Color','k');
    text(first+(T-first)*7/10,q(ti),sprintf('Posterior %3.2f Quantile',tau(ti)));
  end
  line([first,T],[avg,avg],'LineStyle','-','Color','k');
  text(first+(T-first)*3/10,avg,'Posterior Mean');    
end
if (incTrue)
  line([first,T],[rho0,rho0],'LineStyle','-','Color','r')
  text(first+(T-first)/20,rho0,'True Value');
end
title('\rho');
if (sp>nsub)
  print('-depsc2',sprintf('figures/%s%d',suffix,ns));
  ns = ns+1;
  sp = 1;
end
%print('-depsc2',sprintf('figures/theta%s.eps',suffix));
%if (sp>nsub)
xlim([0 T]);

for k=1:numel(all)
  theta = all{k}.theta(first:T);
  avg = mean(theta);
  q = quantile(theta(post),tau);
  if (k==1)
    if (sp==1) 
      figure;
    end
    subplot(sub1,sub2,sp),plot(t,theta,'Color',color{k},'LineStyle', ...
                               lt(k));
    sp=sp+1;
  else 
    line(t,theta,'LineStyle',lt(k),'Color',color{k});
  end
  for ti=1:numel(tau) 
    line([first,T],[q(ti),q(ti)],'LineStyle','--','Color','k');
    text(first+(T-first)*7/10,q(ti),sprintf('Posterior %3.2f Quantile',tau(ti)));
  end
  line([first,T],[avg,avg],'LineStyle','-','Color','k');
  text(first+(T-first)*3/10,avg,'Posterior Mean');    
end
if (incTrue)
  line([first,T],[theta0,theta0],'LineStyle','-','Color','r')
  text(first+(T-first)/20,theta0,'True Value');
end
title('\theta');
%print('-depsc2',sprintf('figures/theta%s.eps',suffix));
%if (sp>nsub)
xlim([0 T]);
if (sp>nsub)
  print('-depsc2',sprintf('figures/%s%d',suffix,ns));
  ns = ns+1;
  sp = 1;
end

for k=1:numel(all)
  shape = all{k}.shape(first:T);
  avg = mean(shape);
  q = quantile(shape(post),tau);
  if (k==1)
    if (sp==1) 
      figure;
    end
    subplot(sub1,sub2,sp),plot(t,shape,'Color',color{k},'LineStyle', ...
                               lt(k));
    sp=sp+1;
  else 
    line(t,shape,'LineStyle',lt(k),'Color',color{k});
  end
  for ti=1:numel(tau) 
    line([first,T],[q(ti),q(ti)],'LineStyle','--','Color','k');
    text(first+(T-first)*7/10,q(ti),sprintf('Posterior %3.2f Quantile',tau(ti)));
  end
  line([first,T],[avg,avg],'LineStyle','-','Color','k');  text(first+(T-first)*3/10,avg,'Posterior Mean');    
end
if (incTrue)
  line([first,T],[shape0,shape0],'LineStyle','-','Color','r')
  text(first+(T-first)/20,shape0,'True Value');
end
title('shape');
%print('-depsc2',sprintf('figures/shape%s.eps',suffix));
%if (sp>nsub)
xlim([0 T]);
print('-depsc2',sprintf('figures/%s%d',suffix,ns));
ns = ns+1;
sp = 1;
%end

