%% creates tables in march 2011 version of paper
config
clear;
path(path,'../ver2');
path(path,'../ver2/counterfactuals');
config;
load(resultFile);
seed = 1337; %sum(100*clock);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
nSim = 10*data.N;

%% compute posterior means
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
if isfield('alpha',chain)
  alpha = mean(chain.alpha(post,:),1)';
else
  alpha = 0;
end

%% creat table of estimates and standard errors
names = {'\omega','\psi','\mu_{\lambda 1}','\mu_{\lambda 2}', ...
         '\mu_{\lambda 3}'};
if (size(data.x{1},2)==13)
  xnames{1} = {'Constant','family','+ spouse','+ child', ...
               'switch 04','switch 05','switch 06','age', ...
               'female','tenure','income','rsFam','rsFam*1{rsFam>5}','2004'};
elseif (size(data.x{1},2)==14)
  xnames{1} = {'Constant','family','+ spouse','+ child', ...
               'switch 04','switch 05','switch 06','age', ...
               'female','tenure','income','rsFamQ2','rsFamQ3','rsFamQ4','2004'};
else
  xnames{1} = {'Constant','family','+ spouse','+ child', ...
               'switch 04','switch 05','switch 06','age', ...
               'female','tenure','income','rsFam','2004'};
end
xnames{2} = xnames{1};
xnames{3} = xnames{1};
options.sampleAll = true;
tic
simdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
toc

%% latex version
out = fopen(sprintf('tex/tables/%sParmPse.tex',prefix),'w');
fprintf(out,'Mean Shifters \\\\ \n');
if sum(abs(alpha))<1e-15
  haveAlpha = false;
  fprintf(out,'& $\\mu_\\lambda$ & $\\kappa_\\lambda$ & $\\omega$ & $\\psi$ \\\\ \n');
  fprintf(out,['& (Health risk)& (Health risk)& (Moral hazard) & (Risk ' ...
               'aversion) \\\\ \n']);
  bindex = [3 0 1 2];
else
  haveAlpha = true;
  fprintf(out,'& $\\mu_\\lambda$ & $\\kappa_\\lambda$ & $\\omega$ & $\\psi$ & $\\alpha$ \\\\ \n');
  fprintf(out,['& (Health risk)& (Health risk)& (Moral hazard) & (Risk ' ...
               'aversion) & (Heteroskedasticity) \\\\ \n']);
  bindex = [3 0 1 2 -1];
end

for j=1:numel(xnames{1})
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(bindex))
    col = j;
    if (bindex(k)==0)  
      vals = -chain.bll(post,:);
    elseif bindex(k)==-1
      col = j-1;
      vals = chain.alpha(post,:);
    else 
      vals = chain.beta{bindex(k)}(post,:);
    end
    if (col<=size(vals,2) && col>0)
      fprintf(out,'& %5.3g (%.3g) ', [mean(vals(:,col)) std(vals(:,col))]);
    end
  end
  fprintf(out,'\\\\ \n');
end  

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Conditional variance-covariance matrix}\\\\ \n');
fprintf(out,'& $\\bar{\\mu}_\\lambda$ & $\\omega$ & $\\psi$\\\\ \n');
fprintf(out,['& (Health risk)& (Moral hazard) & (Risk aversion) \\\\ \n']);
bindex = [3 1 2];
fprintf(out,'$\\bar{\\mu}_\\lambda$ & %.3g (%.3g) & %.3g (%.3g)& %.3g (%.3g)\\\\ \n', ...
        mean(chain.sig(post,bindex(1),4)), ...
        std(chain.sig(post,bindex(1),4)), ...
        mean(chain.sig(post,bindex(1),bindex(2))), ...
        std(chain.sig(post,bindex(1),bindex(2))), ...
        mean(chain.sig(post,bindex(1),bindex(3))), ...
        std(chain.sig(post,bindex(1),bindex(3))));
fprintf(out,'$\\omega$ & -- & %.3g (%.3g)& %.3g (%.3g)\\\\ \n', ...
        mean(chain.sig(post,bindex(2),bindex(2))), ...
        std(chain.sig(post,bindex(2),bindex(2))), ...
        mean(chain.sig(post,bindex(2),bindex(3))), ...
        std(chain.sig(post,bindex(2),bindex(3))));
fprintf(out,'$\\psi$ & -- & -- & %.3g (%.3g)\\\\ \n', ...
        mean(chain.sig(post,bindex(3),bindex(3))), ...
        std(chain.sig(post,bindex(3),bindex(3))));

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Additional parameters}\\\\ \n');
fprintf(out,'$\\sigma_\\epsilon^2$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.sig(post,3,3)-chain.sig(post,3,4)),std(chain.sig(post,3,3)-chain.sig(post,3,4)));
fprintf(out,'$\\sigma_\\mu^2=\\sigma_{\\bar{\\mu}}^2+\\sigma_\\epsilon^2$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.sig(post,3,3)),std(chain.sig(post,3,3)));
fprintf(out,'$\\rho=corr(\\mu_{\\lambda,i1},\\mu_{\\lambda,i2}|X,\\omega,\\psi)$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.rho(post)),std(chain.rho(post)));
fprintf(out,['$corr\\binom{\\mu_{\\lambda,i1}-x_{i1}\\beta_\\lambda,}{\\,\\mu_{\\lambda,i2}-x_{i2}\\beta_\\lambda} = ' ....
             '\\frac{\\sigma_{\\bar{\\mu}}^2}{\\sigma_\\mu^2}$ & %.3g (%.3g)\\\\ \n'], ... 
        mean( (chain.sig(post,3,4)./chain.sig(post,3,3) )), ...
        std( (chain.sig(post,3,4)./chain.sig(post,3,3) )));
fprintf(out,'$\\sigma_\\kappa$& %.3g (%.3g)\\\\ \n', ...
        mean(chain.sigll(post)),std(chain.sigll(post)));
fprintf(out,'$\\gamma_1$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.theta(post)), std(chain.theta(post)));
fprintf(out,'$\\gamma_2$ & %.3g (%.3g)\\\\ \n',...
        mean(chain.shape(post)), std(chain.shape(post)));

%% compute implied quantities and their se's
S = 100; % number of simulations
for s=1:S
  r = ceil(rand()*nSaved/2 + nSaved/2);
  Sigma = squeeze((chain.sig(r,:,:)));
  gamma = squeeze((chain.gamma(r,:,:)));
  theta = squeeze((chain.theta(r)));
  shape = squeeze((chain.shape(r)));
  rho = squeeze((chain.rho(r)));
  bll = (chain.bll(r,:))';
  sigll = (chain.sigll(r));
  if isfield('alpha',chain)
    alpha = (chain.alpha(r,:))';
  else
    alpha = 0;
  end
  sdat{s} = simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                     options);
  meanlatent(s,:) = [mean(sdat{s}.lambda(:)), mean(sdat{s}.omega), ...
                mean(sdat{s}.psi) mean(sdat{s}.eLambda(:))];
  stdlatent(s,:) = [std(sdat{s}.lambda(:)), std(sdat{s}.omega), std(sdat{s}.psi) ...
                   std(sdat{s}.eLambda(:))]; 
  corrlatent(s,:,:) =  [mycorr(sdat{s}.eLambda,sdat{s}.eLambda), ...
                      mycorr(sdat{s}.eLambda,sdat{s}.omega), ...
                      mycorr(sdat{s}.eLambda,sdat{s}.psi); ...
                      NaN, mycorr(sdat{s}.omega,sdat{s}.omega), ...
                      mycorr(sdat{s}.omega,sdat{s}.psi); ...
                      NaN, NaN, mycorr(sdat{s}.psi,sdat{s}.psi)];
  latent = {sdat{s}.eLambda, sdat{s}.omega, sdat{s}.psi };
  for j=1:size(sdat{s}.x{1},2)
    for k=1:(numel(latent))
      if (k==1) 
        margCorr(s,k,j) = mycorr(latent{k}(2,:)',sdat{s}.x{3}{2}(:,j));
      else 
        margCorr(s,k,j) = mycorr(latent{k},sdat{s}.x{1}(:,j));
      end
    end
  end  
  
  y = [sdat{s}.eLambda(1,:)'; sdat{s}.eLambda(2,:)'];
  X = [sdat{s}.x{3}{1}; sdat{s}.x{3}{2}];
  rcoef{1}(s,:) = (X'*X) \ X'*y;
  r2lev(s,1) = var(X*rcoef{1}(s,:)')/var(y);
  for j=2:numel(latent)
    X = sdat{s}.x{j-1};
    rcoef{j}(s,:) = (X'*X) \ X'*latent{j};
    r2lev(s,j) = var(X*rcoef{j}(s,:)')/var(latent{j});
  end
  
  minel = min(sdat{s}.eLambda(:));
  if (minel > 0) 
    minel = 0;
  end
  y = log([sdat{s}.eLambda(1,:)'; sdat{s}.eLambda(2,:)'] + abs(minel) +1);
  X = [sdat{s}.x{3}{1}; sdat{s}.x{3}{2}];
  rcoeflog{1}(s,:) = (X'*X) \ X'*y;
  r2log(s,1) = var(X*rcoeflog{1}(s,:)')/var(y);
  for j=2:numel(latent)
    X = sdat{s}.x{j-1};
    rcoeflog{j}(s,:) = (X'*X) \ X'*log(latent{j});
    r2log(s,j) = var(X*rcoeflog{j}(s,:)')/var(latent{j});
  end

  fprintf('%d of %d finished\n',s,S);
end

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Implied Quantities}\\\\ \n');
fprintf(out,'& $\\lambda$ & $\\omega$ & $\\psi$ & $E[\\lambda]$ \\\\ \n');
fprintf(out,['& (Health risk)& (Moral hazard) & (Risk aversion) & (Health ' ...
             'risk) \\\\ \n']);
fprintf(out,'Expected & %.3g (%.3g) & %.3g  (%.3g) & %.3g  (%.3g) & %.3g  (%.3g) \\\\ \n', ...
        mean(simdata.lambda(:)), std(meanlatent(:,1)), ...
        mean(simdata.omega),  std(meanlatent(:,2)), ...
        mean(simdata.psi),  std(meanlatent(:,3)), ...
        mean(simdata.eLambda(:)),  std(meanlatent(:,4)));
fprintf(out,'Std. Dev. & %.3g (%.3g) & %.3g (%.3g) & %.3g (%.3g) & %.3g (%.3g) \\\\ \n', ...
        std(simdata.lambda(:)), std(stdlatent(:,1)), ...
        std(simdata.omega), std(stdlatent(:,2)), ...
        std(simdata.psi), std(stdlatent(:,3)), ...
        std(simdata.eLambda(:)), std(stdlatent(:,4)));
fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Unconditional correlation matrix} \\\\ \n');
fprintf(out,'& $E[\\lambda]$ & $\\omega$ & $\\psi$ \\\\ \n');
fprintf(out,' $E[\\lambda]$ & %.3g  (%.3g) & %.3g  (%.3g) & %.3g  (%.3g) \\\\ \n', ...
        mycorr(simdata.eLambda,simdata.eLambda), std(corrlatent(:,1,1)), ...
        mycorr(simdata.eLambda,simdata.omega), std(corrlatent(:,1,2)), ...
        mycorr(simdata.eLambda,simdata.psi), std(corrlatent(:,1,3)));
fprintf(out,' $\\omega$ & -- & %.3g (%.3g) & %.3g (%.3g) \\\\ \n', ...
        mycorr(simdata.omega,simdata.omega),std(corrlatent(:,2,2)), ...
        mycorr(simdata.omega,simdata.psi), std(corrlatent(:,2,3)));
fprintf(out,' $\\psi$ & -- & -- & %.3g (%.3g) \\\\ \n', ...
        mycorr(simdata.psi,simdata.psi), std(corrlatent(:,3,3)));

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Marginal correlations} \\\\ \n');
fprintf(out,'& $E[\\lambda]$ & $\\omega$ & $\\psi$ \\\\ \n');
latent = {simdata.eLambda, simdata.omega, simdata.psi };
for j=1:size(simdata.x{1},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    if (k==1) 
      fprintf(out,'& %.3g (%.3g)', mycorr(latent{k}(2,:)',simdata.x{3}{2}(:,j)), ...
              std(margCorr(:,k,j)));
    else 
      fprintf(out,'& %.3g (%.3g)', mycorr(latent{k},simdata.x{1}(:,j)), ...
              std(margCorr(:,k,j)));
    end
  end
  fprintf(out,'\\\\ \n');
end  

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Regression coefficients (levels)} \\\\ \n');
fprintf(out,'& $E[\\lambda]$ & $\\omega$ & $\\psi$ \\\\ \n');
latent = {simdata.eLambda, simdata.omega, simdata.psi};
if (data.T~=2) 
  warning('assuming data.T=2');
end
y = [simdata.eLambda(1,:)'; simdata.eLambda(2,:)'];
X = [simdata.x{3}{1}; simdata.x{3}{2}];
coef{1} = (X'*X) \ X'*y;
r2(1) = var(X*coef{1})/var(y);
for j=2:numel(latent)
  X = simdata.x{j-1};
  coef{j} = (X'*X) \ X'*latent{j};
  r2(j) = var(X*coef{j})/var(latent{j});
end

for j=1:size(simdata.x{3}{2},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    if (j<=numel(coef{k}))
      fprintf(out,'& %.3g (%.3g)', coef{k}(j),std(rcoef{k}(:,j)));
    else 
      fprintf(out,'&      ');
    end
  end
  fprintf(out,'\\\\ \n');
end  
fprintf(out,'$R^2$ & %.3g (%.3g) & %.3g (%.3g) & %.3g (%.3g) \\\\ \n', ...
        r2(1),std(r2lev(:,1)), r2(2),std(r2lev(:,2)), r2(3),std(r2lev(:,3)));

fprintf(out,['\\\\ \n\\multicolumn{4}{l}{Regression coefficients (logs)} ' ...
             '\\\\ \n']);
minel = min(simdata.eLambda(:));
if (minel > 0) 
  minel = 0;
end
fprintf(out,['& $\\log \\begin{pmatrix}E[\\lambda] + 1 \\\\ ' ...
             '    \\vert \\min E[\\lambda]=%.2g\\vert \\end{pmatrix} $' ...
             '& $\\log \\omega$ &  $\\log \\psi$ \\\\ \n'],minel); 
latent = {log(simdata.eLambda + abs(minel)+1), log(simdata.omega), log(simdata.psi)};
if (data.T~=2) 
  warning('assuming data.T=2');
end
y = log([simdata.eLambda(1,:)'; simdata.eLambda(2,:)'] + abs(minel)+1);
X = [simdata.x{3}{1}; simdata.x{3}{2}];
coef{1} = (X'*X) \ X'*y;
r2(1) = var(X*coef{1})/var(y);
for j=2:numel(latent)
  X = simdata.x{j-1};
  coef{j} = (X'*X) \ X'*latent{j};
  r2(j) = var(X*coef{j})/var(latent{j});
end

for j=1:size(simdata.x{3}{2},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    if (j<=numel(coef{k}))
      fprintf(out,'& %.3g (%.3g)', coef{k}(j),std(rcoeflog{k}(:,j)));
    else 
      fprintf(out,'&      ');
    end
  end
  fprintf(out,'\\\\ \n');
end  
fprintf(out,'$R^2$ & %.3g (%.3g) & %.3g (%.3g) & %.3g (%.3g) \\\\ \n', ...
        r2(1),std(r2log(:,1)), r2(2),std(r2log(:,2)), r2(3),std(r2log(:,3)));
fclose(out);


