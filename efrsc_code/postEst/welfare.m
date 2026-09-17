close all;
path(path,'..');
config;
subsidy = 0.0;
prefix = [prefix 's00'];
load(resultFile); 
%prefix = 'moreXrs';
seed = 1337; %sum(100*clock);
setSeed(seed)
if (isOctave)
  debug_on_error(1);
else
  dbstop if error;
end
clear options;
options.cf = true;
options.nosel = true;
options.sampleAll=true;
% compute posterior means
nSaved = numel(chain.meanL);
post = floor(nSaved/2):nSaved;
for k=1:numel(chain.beta);
  beta{k} = mean(chain.beta{k}(post,:),1)';
end
Sigma = squeeze(mean(chain.sig(post,:,:),1));
gamma = squeeze(mean(chain.gamma(post,:,:),1));
theta = squeeze(mean(chain.theta(post)));
shape = squeeze(mean(chain.shape(post)));
rho = squeeze(mean(chain.rho(post)));
bll = mean(chain.bll(post,:),1)';
sigll = mean(chain.sigll(post));
if isfield(chain,'alpha')
  alpha = mean(chain.alpha(post,:),1)';
else
  alpha = 0;
end


% load data without plan 1 deleted
dataEst = data;
negLambdaZeroUtil = false;
if (negLambdaZeroUtil)
  data = loadData('al.csv',true,2004);
  data.order = 20;
  data.varLo = 0;
  data.varHi = 4*var(log(1+data.totalSpend(~isnan(data.totalSpend))));
  avgSpend = mean(data.totalSpend(~isnan(data.totalSpend)));
  data.logOmegaHi = Inf;
  data.lamloHi = 3000;
  data.maxIter = dataEst.maxIter;
  data.order = dataEst.order;
  data.threads = dataEst.threads;
  sel = squeeze(data.avail(5,2,:)==1);
  data = dropBad(data,data.choice(2,:)'==1 & sel);
else 
  dataEst = data;
end
nSim = 10*data.N;

fullTable = [];
[order xint wint] = gqzero(data.order);
xint = xint'*sqrt(2);
order = numel(xint);
wint = wint/sqrt(pi);

setSeed(seed);
options.balance15 = 0;
options.avgPlan = true;
simdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
fprintf('simdata done\n');
[eCost eSpend eOop] = exCost(simdata);
[eCostNMH eSpendNMH] = exCost(simdata,true);

%simdata.exSpend;
%simdata.value;
% how much we need to increase price of 1 by to make each person indifferent
dprem = squeeze(simdata.value(1,2,:)-simdata.value(5,2,:)) + ...
        subsidy*(squeeze(simdata.prem(1,2,:)- simdata.prem(5,2,:)));

% 

out = fopen(sprintf('tex/tables/%sWelfare.tex',prefix),'w');
fprintf(out,['Tier & E[Prem. 5 $-$ Prem. 1] & $P($choose $5)$ & $E[$spend$]$ ' ...
             '&  $E[$total welfare$]$ & $E[\\omega|1]$ & $E[\\omega|5]$ & ' ...
             '$E[\\lambda|1]$ & $E[\\lambda|5]$ \\\\ ' ...
             '\\hline \n']); 
fprintf(out,['\\multicolumn{9}{c}{Asymmetric information baseline} \\\\ ' ...
             '\\hline \n']);
tiername = {'single','family','+spouse','+child'};
%v0 = squeeze(simdata.value(1,2,:) - eCost(1,2,:));
for g=1:4
  weight(g) = sum(simdata.covg==g);
end
weight = weight/sum(weight);

for g=1:4
  dp = findcfPrice(squeeze(eCost([1 5],2,simdata.covg==g)),dprem(simdata.covg==g),subsidy);
  premdiff = unique(simdata.prem(5,2,simdata.covg==g)) - ...
      unique(simdata.prem(1,2,simdata.covg==g)) - dp;
  es = squeeze(eSpend(1,2,:));
  es(dprem<=dp*(1-subsidy)) = squeeze(eSpend(5,2,dprem<=dp*(1-subsidy)));
  ev = squeeze(simdata.value(1,2,:));
  ev(dprem<=dp*(1-subsidy)) = squeeze(simdata.value(5,2,dprem<=dp*(1-subsidy)));
  ev = ev; %- v0;
  ec = squeeze(eCost(1,2,:));
  ec(dprem<=dp*(1-subsidy)) = squeeze(eCost(5,2,dprem<=dp*(1-subsidy)));
  etw = ev - ec;
  
  c5 = dprem<=dp*(1-subsidy);
  table(g,:) = [premdiff, ...
          mean(dprem(simdata.covg==g)<=dp*(1-subsidy)), ...
          mean(es(simdata.covg==g)),mean(etw(simdata.covg==g)), ...
          mean(simdata.omega(simdata.covg==g & dprem>dp*(1-subsidy))), ...
          mean(simdata.omega(simdata.covg==g & dprem<=dp*(1-subsidy))), ...
          mean(simdata.lambda(2,simdata.covg==g & dprem>dp*(1-subsidy))), ...
          mean(simdata.lambda(2,simdata.covg==g & dprem<=dp*(1-subsidy)))];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
          tiername{g}, table(g,:));
  %  for asym info choice, nmh spend
  %premdiff = mean(eCostNMH(5,2,simdata.covg==g & dprem<=dp*(1-subsidy)),3) - ...
  %    mean(eCostNMH(1,2,simdata.covg==g & dprem>dp*(1-subsidy)),3) ...
  %    - unique(simdata.prem(1,2,simdata.covg==g)) + ...
  %    unique(simdata.prem(5,2,simdata.covg==g));
  es = squeeze(eSpendNMH(1,2,:));
  es(dprem<=dp*(1-subsidy)) = squeeze(eSpendNMH(5,2,dprem<=dp*(1-subsidy)));
  ev = squeeze(simdata.valueNMH(1,2,:));
  ev(dprem<=dp*(1-subsidy)) = squeeze(simdata.valueNMH(5,2,dprem<=dp*(1-subsidy)));
  ev = ev; %- v0;
  ec = squeeze(eCostNMH(1,2,:));
  ec(dprem<=dp*(1-subsidy)) = squeeze(eCostNMH(5,2,dprem<=dp*(1-subsidy)));
  etw = ev - ec;
  table2(g,:) = [premdiff, ...
                 mean(dprem(simdata.covg==g)<=dp*(1-subsidy)), ...
                 mean(es(simdata.covg==g)),mean(etw(simdata.covg==g)), ...
                 mean(simdata.omega(simdata.covg==g & dprem>dp*(1-subsidy))), ...
                 mean(simdata.omega(simdata.covg==g & dprem<=dp*(1-subsidy))), ...
                 mean(simdata.lambda(2,simdata.covg==g & dprem>dp*(1-subsidy))), ...
                 mean(simdata.lambda(2,simdata.covg==g & dprem<=dp*(1-subsidy)))];
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
        'All', weight*table);
fullTable = [table; weight*table];
clear table;

fprintf(out,['\\hline \\multicolumn{9}{c}{Asymmetric info choices, no moral hazard spending} \\\\ ' ...
             '\\hline \n']);
for g=1:4
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
          tiername{g}, table2(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
        'All', weight*table2);
fullTable = [fullTable; table2; weight*table2];
clear table2;

%% no moral hazard
fprintf(out,['\\hline \\multicolumn{9}{c}{No moral hazard} \\\\ ' ...
             '\\hline \n']);
dprem = squeeze(simdata.valueNMH(1,2,:)-simdata.valueNMH(5,2,:)) + ...
        subsidy*(squeeze(simdata.prem(1,2,:)- simdata.prem(5,2,:)));

for g=1:4
  dp = findcfPrice(squeeze(eCostNMH([1 5],2,simdata.covg==g)),dprem(simdata.covg==g),subsidy);
  premdiff = unique(simdata.prem(5,2,simdata.covg==g)) - ...
      unique(simdata.prem(1,2,simdata.covg==g)) - dp;
  es = squeeze(eSpendNMH(1,2,:));
  es(dprem<=dp*(1-subsidy)) = squeeze(eSpendNMH(5,2,dprem<=dp*(1-subsidy)));
  ev = squeeze(simdata.valueNMH(1,2,:));
  ev(dprem<=dp*(1-subsidy)) = squeeze(simdata.valueNMH(5,2,dprem<=dp*(1-subsidy)));
  ev = ev; % - v0;
  ec = squeeze(eCostNMH(1,2,:));
  ec(dprem<=dp*(1-subsidy)) = squeeze(eCostNMH(5,2,dprem<=dp*(1-subsidy)));
  etw = ev - ec;
  
  c5 = dprem<=dp*(1-subsidy);
  table(g,:) = [premdiff, ...
                mean(dprem(simdata.covg==g)<=dp*(1-subsidy)), ...
                mean(es(simdata.covg==g)),mean(etw(simdata.covg==g)), ...
                mean(simdata.omega(simdata.covg==g & dprem>dp*(1-subsidy))), ...
                mean(simdata.omega(simdata.covg==g & dprem<=dp*(1-subsidy))), ...
                mean(simdata.lambda(2,simdata.covg==g & dprem>dp*(1-subsidy))), ...
                mean(simdata.lambda(2,simdata.covg==g & dprem<=dp*(1-subsidy)))];          
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
          tiername{g}, table(g,:));

  %  for nmh choice, mh spend
  %premdiff = mean(eCost(5,2,simdata.covg==g & dprem<=dp*(1-subsidy)),3) - ...
  %    mean(eCost(1,2,simdata.covg==g & dprem>dp*(1-subsidy)),3) ...
  %    - unique(simdata.prem(1,2,simdata.covg==g)) + ...
  %    unique(simdata.prem(5,2,simdata.covg==g));
  es = squeeze(eSpend(1,2,:));
  es(dprem<=dp*(1-subsidy)) = squeeze(eSpend(5,2,dprem<=dp*(1-subsidy)));
  ev = squeeze(simdata.value(1,2,:));
  ev(dprem<=dp*(1-subsidy)) = squeeze(simdata.value(5,2,dprem<=dp*(1-subsidy)));
  ev = ev; %- v0;
  ec = squeeze(eCost(1,2,:));
  ec(dprem<=dp*(1-subsidy)) = squeeze(eCost(5,2,dprem<=dp*(1-subsidy)));
  etw = ev - ec;
  table2(g,:) = [premdiff, ...
                 mean(dprem(simdata.covg==g)<=dp*(1-subsidy)), ...
                 mean(es(simdata.covg==g)),mean(etw(simdata.covg==g)), ...
                 mean(simdata.omega(simdata.covg==g & dprem>dp*(1-subsidy))), ...
                 mean(simdata.omega(simdata.covg==g & dprem<=dp*(1-subsidy))), ...
                 mean(simdata.lambda(2,simdata.covg==g & dprem>dp*(1-subsidy))), ...
                 mean(simdata.lambda(2,simdata.covg==g & dprem<=dp*(1-subsidy)))];
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];
clear table;

fprintf(out,['\\hline \\multicolumn{9}{c}{No moral hazard choices, spending with moral hazard} \\\\ ' ...
             '\\hline \n']);
for g=1:4
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
          tiername{g}, table2(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g \\\\ \n'], ...
        'All', weight*table2);
fullTable = [fullTable; table2; weight*table2];
clear table2;


fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about $E[\\lambda]$} \\\\ ' ...
             '\\hline \n']);
%% welfare with symmetric information about lambda
% repeat same as above, but find kernel weighted zero for each
% observation
dprem = (squeeze(simdata.value(1,2,:)-simdata.value(5,2,:))) + ...
        subsidy*(squeeze(simdata.prem(1,2,:)- simdata.prem(5,2,:)));

x = log(simdata.eSpendNI+1);
for g=1:4
  ind = find(simdata.covg==g);
  %ind = ind(randperm(numel(ind)));
  N = numel(ind);
  es = zeros(N,1);
  ev = zeros(N,1);
  ec = zeros(N,1);
  ch = zeros(N,1);
  etw = zeros(N,1);
  %xi = quantile(x(ind),0.025:0.025:0.975);
  dq = 0.02;
  %tau = dq:dq:(1-dq);
  tau = 0:dq:(1);
  xi = unique(quantile(x(ind),tau));
  if (numel(xi) == 1)
    xi = [xi - 1, xi + 1];
  end
  [junk rankOrder] = sort(x(ind));
  tauL = tau-0.1;
  tauH = tau+0.1;
  qx = rankOrder/numel(x(ind));  
  clear dpi dwi dvi dci;
  for i=1:numel(xi)
    % Silverman's rule for bandwidth
    bw = 0.9 * std(x(simdata.covg==g))*sum(simdata.covg==g)^(-1/5); 
    bw = 3*bw;
    w = exp(-.5*((x(ind) - xi(i))/bw).^2)/ ...
        sqrt(2*pi);
    %w = qx>tauL(i) & qx<tauH(i);
    dpi(i) = findcfPrice(squeeze(eCost([1 5],2,ind)),dprem(ind),subsidy,w);

    ecw = squeeze(eCost([1 5], 2, ind)).*(ones(2,1)*w');
    evw = squeeze(simdata.value([1 5], 2, ind)).*(ones(2,1)*w');
    eww = evw - ecw;        
    dwi(i) = mean(eww(2,:)-eww(1,:))/mean(w);
    dvi(i) = mean(evw(2,:)-evw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
    dci(i) = mean(ecw(2,:)-ecw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
  end
  % linearly interpolate
  dp = interp1(xi,dpi,x,'pchip','extrap');
  premdiff = unique(simdata.prem(5,2,simdata.covg==g)) - ...
      unique(simdata.prem(1,2,simdata.covg==g)) - dp;
  es = squeeze(eSpend(1,2,ind));
  ev = squeeze(simdata.value(1,2,ind));
  ec = squeeze(eCost(1,2,ind));
  ch(:)  = 1;
  ch(dp(ind)*(1-subsidy)>dprem(ind)) = 5;
  es(ch==5) = squeeze(eSpend(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(simdata.value(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCost(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(premdiff(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
  figure;
  if (numel(ind)>300)
    ind = ind(1:floor(numel(ind)/300):end);
  end
  dv = interp1(xi,dvi,x,'pchip','extrap');
  dc = interp1(xi,dci,x,'pchip','extrap');
  dw = interp1(xi,dwi,x,'pchip','extrap');
  plot(x(ind),premdiff(ind),'.',x(ind),dv(ind),'.', ...
       x(ind),dc(ind),'.',x(ind),dw(ind),'.');
  xlabel('log(E[\lambda])');
  legend('premium difference','willingness to pay', ...
         'cost difference','E[Welfare(5)-Welfare(1)]');
  %print(sprintf('figures/%sMuLPrem%d.tiff',prefix,g),'-dtiff');
  print(sprintf('figures/%sMuLPrem%d.eps',prefix,g),'-depsc2');
  %system(sprintf('tiff2ps -e -3 figures/%sMuLPrem%d.tiff > figures/%sMuLPrem%d.eps', ...
  %               prefix,g,prefix,g));
  fprintf('%d done\n',g);
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about $E[\\lambda]$ and no moral hazard}\\\\ ' ...
             '\\hline \n']);
%% welfare with symmetric information about lambda and no moral hazard
% repeat same as above, but find kernel weighted zero for each
% observation
dprem = (squeeze(simdata.valueNMH(1,2,:)-simdata.valueNMH(5,2,:))) + ...
        subsidy*(squeeze(simdata.prem(1,2,:)- simdata.prem(5,2,:)));

x = log(simdata.eSpendNI+1);
for g=1:4
  ind = find(simdata.covg==g);
  %ind = ind(randperm(numel(ind)));
  N = numel(ind);
  es = zeros(N,1);
  ev = zeros(N,1);
  ec = zeros(N,1);
  ch = zeros(N,1);
  etw = zeros(N,1);
  %xi = quantile(x(ind),0.025:0.025:0.975);
  dq = 0.02;
  %tau = dq:dq:(1-dq);
  tau = 0:dq:(1);
  xi = unique(quantile(x(ind),tau));
  if (numel(xi) == 1)
    xi = [xi - 1, xi + 1];
  end
  [junk rankOrder] = sort(x(ind));
  tauL = tau-0.1;
  tauH = tau+0.1;
  qx = rankOrder/numel(x(ind));  
  clear dpi dwi dvi dci;
  for i=1:numel(xi)
        % Silverman's rule for bandwidth 
    bw = 0.9 * std(x(simdata.covg==g))*sum(simdata.covg==g)^(-1/5); 
    bw =3*bw;
    w = exp(-.5*((x(ind) - xi(i))/bw).^2)/ ...
        sqrt(2*pi);
    %w = qx>tauL(i) & qx<tauH(i);
    dpi(i) = findcfPrice(squeeze(eCostNMH([1 5],2,ind)),dprem(ind), ...
                         subsidy,w);
    
    ecw = squeeze(eCostNMH([1 5], 2, ind)).*(ones(2,1)*w');
    evw = squeeze(simdata.valueNMH([1 5], 2, ind)).*(ones(2,1)*w');
    eww = evw - ecw;    
    dwi(i) = mean(eww(2,:)-eww(1,:))/mean(w);
    dvi(i) = mean(evw(2,:)-evw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
    dci(i) = mean(ecw(2,:)-ecw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
  end
  % linearly interpolate
  dp = interp1(xi,dpi,x,'pchip','extrap');
  premdiff = unique(simdata.prem(5,2,simdata.covg==g)) - ...
      unique(simdata.prem(1,2,simdata.covg==g)) - dp/(1-subsidy);
  es = squeeze(eSpendNMH(1,2,ind));
  ev = squeeze(simdata.valueNMH(1,2,ind));
  ec = squeeze(eCostNMH(1,2,ind));
  ch(:)  = 1;
  ch(dp(ind)*(1-subsidy)>dprem(ind)) = 5;
  es(ch==5) = squeeze(eSpendNMH(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(simdata.valueNMH(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCostNMH(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(premdiff(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
  figure;
  if (numel(ind)>1000)
    ind = ind(1:floor(numel(ind)/1000):end);
  end
  dv = interp1(xi,dvi,x,'pchip','extrap');
  dc = interp1(xi,dci,x,'pchip','extrap');
  dw = interp1(xi,dwi,x,'pchip','extrap');
  plot(x(ind),premdiff(ind),'.',x(ind),dv(ind),'.', ...
       x(ind),dc(ind),'.',x(ind),dw(ind),'.');
  xlabel('log(E[\lambda])');
  legend('premium difference','willingness to pay', ...
         'cost difference','E[Welfare(5)-Welfare(1)]');
  %print(sprintf('figures/%sMuLPremNMH%d.tiff',prefix,g),'-dtiff');
  print(sprintf('figures/%sMuLPremNMH%d.eps',prefix,g),'-depsc2');
  %system(sprintf('tiff2ps -e -3 figures/%sMuLPremNMH%d.tiff > figures/%sMuLPremNMH%d.eps', ...
  %               prefix,g,prefix,g));

  fprintf('%d done\n',g);
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

%% welfare with symmetric information about mu,kappa,& sigma and no moral hazard
% no asym info at all, so each person faces actuarially fair price
fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about '...
             '$(\\mu_\\lambda,\\sigma_\\lambda,\\kappa)$ and no moral hazard}\\\\ ' ...
             '\\hline \n']);
dprem = squeeze(simdata.valueNMH(5,2,:)-simdata.valueNMH(1,2,:) ...
                + simdata.prem(5,2,:) - simdata.prem(1,2,:));
dp = squeeze(eCostNMH(5,2,:) - eCostNMH(1,2,:) - simdata.prem(1,2,:) + simdata.prem(5,2,:));
ch5 = 1 + 4*(dprem>=dp);
for g=1:4
  ind = find(simdata.covg==g);
  es = squeeze(eSpendNMH(1,2,ind));
  ev = squeeze(simdata.valueNMH(1,2,ind));
  ec = squeeze(eCostNMH(1,2,ind));
  ch = ch5(ind);
  es(ch==5) = squeeze(eSpendNMH(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(simdata.valueNMH(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCostNMH(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(dp(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

c1 = find(ch5==1);
names = {'omega','psi','kappa','sigmaLambda','muLambda','deduct1','maxoop1', ...
         'deduct5','maxoop5','E[cost5-cost1]','V5-V1'};
n = 5;
writecsv(sprintf('csv/%schoose1.csv',prefix),names, ...
         [simdata.omega(c1(1:n)), simdata.psi(c1(1:n)), simdata.lamlo(c1(1:n)), ...
          simdata.sigL(c1(1:n)), simdata.muL(c1(1:n),2), ...
          squeeze(simdata.deduct(1,2,c1(1:n))), ...
          squeeze(simdata.maxoop(1,2,c1(1:n))), ...
          squeeze(simdata.deduct(5,2,c1(1:n))), ...
          squeeze(simdata.maxoop(5,2,c1(1:n))), ...
          dp(c1(1:n)), dprem(c1(1:n)) ]);
fprintf('wrote choose1.csv\n');
%% welfare with symmetric information about mu,kappa,& sigma and no moral hazard
% no asym info at all, so each person faces actuarially fair price
fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about ' ...
             '$(\\mu_\\lambda,\\sigma_\\lambda,\\kappa,\\omega)$ with moral hazard }\\\\ ' ...
             '\\hline \n']);
dprem = squeeze(simdata.value(5,2,:)-simdata.value(1,2,:) ...
                + simdata.prem(5,2,:) - simdata.prem(1,2,:));
dp = squeeze(eCost(5,2,:) - eCost(1,2,:) - simdata.prem(1,2,:) + simdata.prem(5,2,:));
ch5 = 1+4*(dprem>=dp);
for g=1:4
  ind = find(simdata.covg==g);
  es = squeeze(eSpend(1,2,ind));
  ev = squeeze(simdata.value(1,2,ind));
  ec = squeeze(eCost(1,2,ind));
  ch = ch5(ind);
  es(ch==5) = squeeze(eSpend(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(simdata.value(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCost(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(dp(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

fprintf('findingChoices ...');
[c0 value0] =  ...
    findChoices_c(log(simdata.omega),1, simdata.lamlo, 1, simdata.muL', ...
                  simdata.sigL, simdata.psi, ... 
                  xint,wint, simdata.deduct, ...
                  simdata.maxoop, simdata.prem, simdata.avail, ...
                  simdata.choice, simdata.totalSpend, ...
                  simdata.maxIter,0.1,1);
fprintf(' done\n findingChoices ...');
[cNMH0 valueNMH0] =  ...
    findChoicesNMH_c(log(simdata.omega),1, simdata.lamlo, 1, simdata.muL', ...
                     simdata.sigL, simdata.psi, ... 
                     xint,wint, simdata.deduct, ...
                     simdata.maxoop, simdata.prem, simdata.avail, ...
                     simdata.choice, simdata.totalSpend, ...
                     simdata.maxIter,1);
fprintf(' done\n ');
for c=1:5
  value0(c,:,:) = squeeze(value0(c,:,:))./ ...
      (ones(simdata.T,1)*simdata.psi');
  value0(~(simdata.avail==1)) = NaN;
  valueNMH0(c,:,:) = squeeze(valueNMH0(c,:,:))./ ...
      (ones(simdata.T,1)*simdata.psi');
  valueNMH0(~(simdata.avail==1)) = NaN;
end
[eCost eSpend eOop] = exCost(simdata,false,true);
[eCostNMH eSpendNMH] = exCost(simdata,true,true);


%% welfare with symmetric information about mu,kappa,& sigma and no moral hazard
% no asym info at all, so each person faces actuarially fair price
fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about '...
             '$(\\mu_\\lambda,\\sigma_\\lambda,\\kappa)$ and no moral hazard and no negative $\\lambda$}\\\\ ' ...
             '\\hline \n']);
dprem = squeeze(valueNMH0(5,2,:)-valueNMH0(1,2,:) ...
                + simdata.prem(5,2,:) - simdata.prem(1,2,:));
dp = squeeze(eCostNMH(5,2,:) - eCostNMH(1,2,:) - simdata.prem(1,2,:) + simdata.prem(5,2,:));
ch5 = 1 + 4*(dprem>=dp);
for g=1:4
  ind = find(simdata.covg==g);
  es = squeeze(eSpendNMH(1,2,ind));
  ev = squeeze(valueNMH0(1,2,ind));
  ec = squeeze(eCostNMH(1,2,ind));
  ch = ch5(ind);
  es(ch==5) = squeeze(eSpendNMH(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(valueNMH0(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCostNMH(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(dp(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

%% welfare with symmetric information about mu,kappa,& sigma and no moral hazard
% no asym info at all, so each person faces actuarially fair price
fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about ' ...
             '$(\\mu_\\lambda,\\sigma_\\lambda,\\kappa,\\omega)$ with moral hazard and no negative $\\lambda$ }\\\\ ' ...
             '\\hline \n']);
dprem = squeeze(value0(5,2,:)-value0(1,2,:) ...
                + simdata.prem(5,2,:) - simdata.prem(1,2,:));
dp = squeeze(eCost(5,2,:) - eCost(1,2,:) - simdata.prem(1,2,:) + simdata.prem(5,2,:));
ch5 = 1+4*(dprem>=dp);
for g=1:4
  ind = find(simdata.covg==g);
  es = squeeze(eSpend(1,2,ind));
  ev = squeeze(value0(1,2,ind));
  ec = squeeze(eCost(1,2,ind));
  ch = ch5(ind);
  es(ch==5) = squeeze(eSpend(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(value0(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCost(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(dp(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];


%% welfare with symmetric information about omega
% repeat same as above, but find kernel weighted zero for each
% observation
fprintf(out,['\\hline \\multicolumn{9}{c}{Symmetric information about ' ...
             '$(\\omega)$ with moral hazard }\\\\ ' ...
             '\\hline \n']);
dprem = squeeze(simdata.value(1,2,:)-simdata.value(5,2,:)) + ...
        subsidy*(squeeze(simdata.prem(1,2,:)- simdata.prem(5,2,:)));
x = log(simdata.omega);
for g=1:4
  ind = find(simdata.covg==g);
  %ind = ind(randperm(numel(ind)));
  N = numel(ind);
  es = zeros(N,1);
  ev = zeros(N,1);
  ec = zeros(N,1);
  ch = zeros(N,1);
  etw = zeros(N,1);
  %xi = quantile(x(ind),0.025:0.025:0.975);
  dq = 0.02;
  %tau = dq:dq:(1-dq);
  tau = 0:dq:(1);
  xi = unique(quantile(x(ind),tau));
  if (numel(xi) == 1)
    xi = [xi - 1, xi + 1];
  end
  [junk rankOrder] = sort(x(ind));
  tauL = tau-0.1;
  tauH = tau+0.1;
  qx = rankOrder/numel(x(ind));  
  clear dpi dwi dvi dci;
  for i=1:numel(xi)
        % Silverman's rule for bandwidth
    bw = 0.9 * std(x(simdata.covg==g))*sum(simdata.covg==g)^(-1/5); 
    bw = 3*bw;
    w = exp(-.5*((x(ind) - xi(i))/bw).^2)/ ...
        sqrt(2*pi);
    %w = qx>tauL(i) & qx<tauH(i);
    dpi(i) = findcfPrice(squeeze(eCost([1 5],2,ind)),dprem(ind),subsidy,w);
    
    ecw = squeeze(eCost([1 5], 2, ind)).*(ones(2,1)*w');
    evw = squeeze(simdata.value([1 5], 2, ind)).*(ones(2,1)*w');
    eww = evw - ecw;    
    dwi(i) = mean(eww(2,:)-eww(1,:))/mean(w);
    dvi(i) = mean(evw(2,:)-evw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
    dci(i) = mean(ecw(2,:)-ecw(1,:))/mean(w) - ...
             unique(simdata.prem(1,2,simdata.covg==g)) + ...
             unique(simdata.prem(5,2,simdata.covg==g));
  end
  % linearly interpolate
  dp = interp1(xi,dpi,x,'pchip','extrap');
  premdiff = unique(simdata.prem(5,2,simdata.covg==g)) - ...
      unique(simdata.prem(1,2,simdata.covg==g)) - dp;
  es = squeeze(eSpend(1,2,ind));
  ev = squeeze(simdata.value(1,2,ind));
  ec = squeeze(eCost(1,2,ind));
  ch(:)  = 1;
  ch(dp(ind)*(1-subsidy)>dprem(ind)) = 5;
  es(ch==5) = squeeze(eSpend(5,2,ind(ch==5)));
  ev(ch==5) = squeeze(simdata.value(5,2,ind(ch==5)));
  ec(ch==5) = squeeze(eCost(5,2,ind(ch==5)));
  etw = ev - ec;
  table(g,:) = [mean(premdiff(ind)), ...
                mean(ch==5),mean(es),mean(etw), ...
                mean(simdata.omega(ind(ch==1))), ...
                mean(simdata.omega(ind(ch==5))), ...
                mean(simdata.lambda(2,ind(ch==1))'), ...
                mean(simdata.lambda(2,ind(ch==5))')];
  fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
          tiername{g}, table(g,:));
  figure;
  if (numel(ind)>300)
    ind = ind(1:floor(numel(ind)/300):end);
  end
  dv = interp1(xi,dvi,x,'pchip','extrap');
  dc = interp1(xi,dci,x,'pchip','extrap');
  dw = interp1(xi,dwi,x,'pchip','extrap');
  plot(x(ind),premdiff(ind),'.',x(ind),dv(ind),'.', ...
       x(ind),dc(ind),'.',x(ind),dw(ind),'.');
  xlabel('log(E[\lambda])');
  legend('premium difference','willingness to pay', ...
         'cost difference','E[Welfare(5)-Welfare(1)]');
  %print(sprintf('figures/%sMuLPrem%d.tiff',prefix,g),'-dtiff');
  print(sprintf('figures/%sOmPrem%d.eps',prefix,g),'-depsc2');
  %system(sprintf('tiff2ps -e -3 figures/%sMuLPrem%d.tiff > figures/%sMuLPrem%d.eps', ...
  %               prefix,g,prefix,g));
  fprintf('%d done\n',g);
end
fprintf(out,['%s & %.0f & %.3f & %.0f & %.0f & %.3g & %.3g & %.3g & %.3g\\\\ \n'], ...
        'All', weight*table);
fullTable = [fullTable; table; weight*table];

fprintf(out,['\\hline \\multicolumn{9}{l}{With moral hazard, welfare with '...
             'all select 1 = %.0f, with all select 5 = %.0f}\\\\ ' ...
             '\n'],mean(simdata.value(1,2,:)-eCost(1,2,:)), ...
        mean(simdata.value(5,2,:)-eCost(5,2,:)));
fprintf(out,['\\multicolumn{9}{l}{Without moral hazard, welfare with '...
             'all select 1 = %.0f, with all select 5 = %.0f}\\\\ ' ...
             '\\hline \n'],mean(simdata.valueNMH(1,2,:)-eCostNMH(1,2,:)), ...
        mean(simdata.valueNMH(5,2,:)-eCostNMH(5,2,:)));

fclose(out);

save(sprintf('csv/%sWelfareTable',prefix),'fullTable');

%% create table that matches paper version
out = fopen(sprintf('tex/tables/%sWelfareP.tex',prefix),'w');
fprintf(out,[' & & Average equilibrium  (incremental) premium & No deductible  ' ...
             'share & Expected spending  per employee & Total welfare  per ' ...
             'employee \\\\ \n']);
fprintf(out,' (1) & "Status quo":  no screening or monitoring ');
r = 5;
fprintf(out,'& %5.0f & %.2f & %5.0f & normalized to 0 \\\\ \n', ...
        fullTable(r,1), fullTable(r,2), fullTable(r,3) );
fprintf(out,[' (2) & "Perfect screening":  premiums depend on $F(\\lambda)$ and ' ...
             '$\\omega$']);
r = 40;
fprintf(out,'& %5.0f & %.2f & %5.0f & %5.0f \\\\ \n', ...
        fullTable(r,1), fullTable(r,2), fullTable(r,3), fullTable(r,4)-fullTable(5,4) ); 
fprintf(out,[' (3) & "Imperfect screening":  premiums depend on $\\omega$ (' ...
             'not $F(\\lambda)$)' ]);
r = 55;
fprintf(out,'& %5.0f & %.2f & %5.0f & %5.0f \\\\ \n', ...
        fullTable(r,1), fullTable(r,2), fullTable(r,3), fullTable(r,4)-fullTable(5,4) ); 
fprintf(out,[' (4) & "Perfect monitoring": contracts only reimburse ' ...
             '$\\lambda$']);
r = 15;
fprintf(out,'& %5.0f & %.2f & %5.0f & %5.0f \\\\ \n', ...
        fullTable(r,1), fullTable(r,2), fullTable(r,3), fullTable(r,4)- fullTable(5,4) ); 
r = 20;
fprintf(out,['(5) & "Imperfect monitoring": perfect monitoring for choice, ' ...
             'but not utilization']);
fprintf(out,'& %5.0f & %.2f & %5.0f & %5.0f \\\\ \n', ...
        fullTable(r,1), fullTable(r,2), fullTable(r,3), fullTable(r,4)-fullTable(5,4) ); 

fclose(out);
