%% creates tables in march 2011 version of paper
clear;
path(path,'..');
config;
load(resultFile);
seed = 1337; %sum(100*clock);
setSeed(seed);
nSim = data.N*10;

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
if isfield(chain,'alpha')
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
options.multMH = false;
tic
simdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
toc

%% table of estimates
out = fopen(sprintf('csv/%sParm.csv',prefix),'w');
fprintf(out,'Mean Shifters, \n');
fprintf(out,', Mu_Lambda , Kappa_Lambda , Omega , Psi\n');
fprintf(out,[', (Health risk), (Health risk), (Moral hazard) , (Risk ' ...
             'aversion) \n']);
bindex = [3 0 1 2];
for j=1:numel(xnames{1})
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(bindex))
    if (bindex(k)==0)  
      vals = chain.bll(post,:);
    else 
      vals = chain.beta{bindex(k)}(post,:);
    end
    if (j<=size(vals,2))
      fprintf(out,', %6.3g (%.3g) ', [mean(vals(:,j)) std(vals(:,j))]);
    end
  end
  fprintf(out,'\n');
end  
fprintf(out,'\nVariance-covariance matrix\n');
fprintf(out,', Mu_Lambda , Omega , Psi\n');
fprintf(out,[', (Health risk), (Moral hazard) , (Risk aversion) \n']);
bindex = [3 1 2];
fprintf(out,'Mu_Lambda , %.3g (%.3g) , %.3g (%.3g), %.3g (%.3g)\n', ...
        mean(chain.sig(post,bindex(1),bindex(1))), ...
        std(chain.sig(post,bindex(1),bindex(1))), ...
        mean(chain.sig(post,bindex(1),bindex(2))), ...
        std(chain.sig(post,bindex(1),bindex(2))), ...
        mean(chain.sig(post,bindex(1),bindex(3))), ...
        std(chain.sig(post,bindex(1),bindex(3))));
fprintf(out,'Omega , -- , %.3g (%.3g), %.3g (%.3g)\n', ...
        mean(chain.sig(post,bindex(2),bindex(2))), ...
        std(chain.sig(post,bindex(2),bindex(2))), ...
        mean(chain.sig(post,bindex(2),bindex(3))), ...
        std(chain.sig(post,bindex(2),bindex(3))));
fprintf(out,'Psi , -- , -- , %.3g (%.3g)\n', ...
        mean(chain.sig(post,bindex(3),bindex(3))), ...
        std(chain.sig(post,bindex(3),bindex(3))));

fprintf(out,'\nAdditional parameters\n');
fprintf(out,'Sigma_Mu_Bar (Rho?), %.3g (%.3g)\n', ...
        mean(chain.rho(post)),std(chain.rho(post)));
fprintf(out,'Sigma_Kappa, %.3g (%.3g)\n', ...
        mean(chain.sigll(post)),std(chain.sigll(post)));
fprintf(out,'Gamma1, %.3g (%.3g)\n', ...
        mean(chain.theta(post)), std(chain.theta(post)));
fprintf(out,'Gamma2, %.3g (%.2g)\n',...
        mean(chain.shape(post)), std(chain.shape(post)));

fprintf(out,'\nImplied Quantities\n');
fprintf(out,', Lambda , Omega , Psi\n');
fprintf(out,', (Health risk), (Moral hazard) , (Risk aversion) \n');
fprintf(out,'Expected , %.3g, %.3g, %.3g \n', ...
        mean(simdata.lambda(:)), mean(simdata.omega), mean(simdata.psi));
fprintf(out,'Std. Dev. , %.3g, %.3g, %.3g \n', ...
        std(simdata.lambda(:)), std(simdata.omega), std(simdata.psi));
fclose(out);


%% latex version
out = fopen(sprintf('tex/tables/%sParmP.tex',prefix),'w');
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
      vals = chain.bll(post,:);
    elseif bindex(k)==-1
      col = j-1;
      vals = chain.alpha(post,:);
    else 
      vals = chain.beta{bindex(k)}(post,:);
    end
    if (col<=size(vals,2) && col>0)
      fprintf(out,'& %5.2g (%.2g) ', [mean(vals(:,col)) std(vals(:,col))]);
    else
      fprintf(out,'& ');
    end
  end
  fprintf(out,'\\\\ \n');
end  
fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Conditional variance-covariance matrix}\\\\ \n');
fprintf(out,'& $\\bar{\\mu}_\\lambda$ & $\\omega$ & $\\psi$\\\\ \n');
fprintf(out,['& (Health risk)& (Moral hazard) & (Risk aversion) \\\\ \n']);
bindex = [3 1 2];
fprintf(out,'$\\bar{\\mu}_\\lambda$ & %.2g (%.2g) & %.2g (%.2g)& %.2g (%.2g)\\\\ \n', ...
        mean(chain.sig(post,bindex(1),4)), ...
        std(chain.sig(post,bindex(1),4)), ...
        mean(chain.sig(post,bindex(1),bindex(2))), ...
        std(chain.sig(post,bindex(1),bindex(2))), ...
        mean(chain.sig(post,bindex(1),bindex(3))), ...
        std(chain.sig(post,bindex(1),bindex(3))));
fprintf(out,'$\\omega$ & -- & %.2g (%.2g)& %.2g (%.2g)\\\\ \n', ...
        mean(chain.sig(post,bindex(2),bindex(2))), ...
        std(chain.sig(post,bindex(2),bindex(2))), ...
        mean(chain.sig(post,bindex(2),bindex(3))), ...
        std(chain.sig(post,bindex(2),bindex(3))));
fprintf(out,'$\\psi$ & -- & -- & %.2g (%.2g)\\\\ \n', ...
        mean(chain.sig(post,bindex(3),bindex(3))), ...
        std(chain.sig(post,bindex(3),bindex(3))));

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Additional parameters}\\\\ \n');
fprintf(out,'$\\sigma_\\epsilon^2$ & %.2g (%.2g)\\\\ \n', ...
        mean(chain.sig(post,3,3)-chain.sig(post,3,4)),std(chain.sig(post,3,3)-chain.sig(post,3,4)));
fprintf(out,'$\\sigma_\\mu^2=\\sigma_{\\bar{\\mu}}^2+\\sigma_\\epsilon^2$ & %.2g (%.2g)\\\\ \n', ...
        mean(chain.sig(post,3,3)),std(chain.sig(post,3,3)));
fprintf(out,'$\\rho=corr(\\mu_{\\lambda,i1},\\mu_{\\lambda,i2}|X,\\omega,\\psi)$ & %.2g (%.2g)\\\\ \n', ...
        mean(chain.rho(post)),std(chain.rho(post)));
fprintf(out,['$corr\\binom{\\mu_{\\lambda,i1}-x_{i1}\\beta_\\lambda,}' ...
             '{\\,\\mu_{\\lambda,i2}-x_{i2}\\beta_\\lambda} = ' ...
             '\\frac{\\sigma_{\\bar{\\mu}}^2}{\\sigma_\\mu^2}$ & %.2g (%.2g)\\\\ \n'], ... 
        mean( (chain.sig(post,3,4)./chain.sig(post,3,3) )), ...
        std( (chain.sig(post,3,4)./chain.sig(post,3,3) )));
fprintf(out,'$\\sigma_\\epsilon$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.rho(post)),std(chain.rho(post)));

fprintf(out,'$\\sigma_kappa$& %.3g (%.3g)\\\\ \n', ...
        mean(chain.sigll(post)),std(chain.sigll(post)));
fprintf(out,'$\\gamma_1$ & %.3g (%.3g)\\\\ \n', ...
        mean(chain.theta(post)), std(chain.theta(post)));
fprintf(out,'$\\gamma_2$ & %.3g (%.3g)\\\\ \n',...
        mean(chain.shape(post)), std(chain.shape(post)));

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Implied Quantities}\\\\ \n');
fprintf(out,'& Lambda & $\\omega$ & $\\psi$\\\\ \n');
fprintf(out,'& (Health risk)& (Moral hazard) & (Risk aversion) \\\\ \n');
fprintf(out,'Expected & %.3g& %.3g& %.3g \\\\ \n', ...
        mean(simdata.lambda(:)), mean(simdata.omega), mean(simdata.psi));
fprintf(out,'Std. Dev. & %.3g& %.3g& %.3g \\\\ \n', ...
        std(simdata.lambda(:)), std(simdata.omega), std(simdata.psi));
fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Unconditional correlation matrix} \\\\ \n');
fprintf(out,'& E[Lambda] & $\\omega$ & $\\psi$\\\\ \n');
fprintf(out,' Lambda & %.3g & %.3g & %.3g \\\\ \n', ...
        mycorr(simdata.eLambda,simdata.eLambda), ...
        mycorr(simdata.eLambda,simdata.omega), ...
        mycorr(simdata.eLambda,simdata.psi));
fprintf(out,' $\\omega$ & -- & %.3g & %.3g \\\\ \n', ...
        mycorr(simdata.omega,simdata.omega), ...
        mycorr(simdata.omega,simdata.psi));
fprintf(out,' $\\psi$ & -- & -- & %.3g \\\\ \n', ...
        mycorr(simdata.psi,simdata.psi));

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Marginal correlations} \\\\ \n');
fprintf(out,'& E[Lambda] & $\\omega$ & $\\psi$ \\\\ \n');
latent = {simdata.eLambda, simdata.omega, simdata.psi };
for j=1:size(simdata.x{1},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    fprintf(out,'& %.2g ', mycorr(latent{k},simdata.x{1}(:,j)));
  end
  fprintf(out,'\\\\ \n');
end  

fprintf(out,'\\\\ \n\\multicolumn{4}{l}{Regression coefficients} \\\\ \n');
fprintf(out,'& E[Lambda] & $\\omega$ & $\\psi$ \\\\ \n');
latent = {simdata.eLambda, simdata.omega, simdata.psi};
if (data.T~=2) 
  warning('assuming data.T=2');
end
y = [simdata.eLambda(1,:)'; simdata.eLambda(2,:)'];
X = [simdata.x{3}{1}; simdata.x{3}{2}];
coef{1} = (X'*X) \ X'*y;
for j=2:numel(latent)
  X = simdata.x{j-1};
  coef{j} = (X'*X) \ X'*latent{j};
end

for j=1:size(simdata.x{3}{2},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    if (j<=numel(coef{k}))
      fprintf(out,'& %.2g ', coef{k}(j));
    else 
      fprintf(out,'&      ');
    end
  end
  fprintf(out,'\\\\ \n');
end  
fclose(out);

