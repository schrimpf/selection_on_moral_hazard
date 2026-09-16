%% csv version
out = fopen(sprintf('csv/%sParmPse.csv',prefix),'w');
fprintf(out,'Mean Shifters  \n');
if sum(abs(alpha))<1e-15
  haveAlpha = false;
  fprintf(out,', MuLambda , KappaLambda , Omega , Psi  \n');
  fprintf(out,[', (Health risk), (Health risk), (Moral hazard) , (Risk ' ...
               'aversion)  \n']);
  bindex = [3 0 1 2];
else
  haveAlpha = true;
  fprintf(out,', MuLambda , KappaLambda , Omega , Psi , Alpha  \n');
  fprintf(out,[', (Health risk), (Health risk), (Moral hazard) , (Risk ' ...
               'aversion) , (Heteroskedasticity)  \n']);
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
      fprintf(out,', %.4g (%.4g) ', [mean(vals(:,col)) std(vals(:,col))]);
    end
  end
  fprintf(out,' \n');
end  
fprintf(out,' \nConditional variance-covariance matrix \n');
fprintf(out,', MuLambda , Omega , Psi \n');
fprintf(out,[', (Health risk), (Moral hazard) , (Risk aversion)  \n']);
bindex = [3 1 2];
fprintf(out,'MuLambda , %.4g (%.4g) , %.4g (%.4g), %.4g (%.4g) \n', ...
        mean(chain.sig(post,bindex(1),4)), ...
        std(chain.sig(post,bindex(1),4)), ...
        mean(chain.sig(post,bindex(1),bindex(2))), ...
        std(chain.sig(post,bindex(1),bindex(2))), ...
        mean(chain.sig(post,bindex(1),bindex(3))), ...
        std(chain.sig(post,bindex(1),bindex(3))));
fprintf(out,'Omega , -- , %.4g (%.4g), %.4g (%.4g) \n', ...
        mean(chain.sig(post,bindex(2),bindex(2))), ...
        std(chain.sig(post,bindex(2),bindex(2))), ...
        mean(chain.sig(post,bindex(2),bindex(3))), ...
        std(chain.sig(post,bindex(2),bindex(3))));
fprintf(out,'Psi , -- , -- , %.4g (%.4g) \n', ...
        mean(chain.sig(post,bindex(3),bindex(3))), ...
        std(chain.sig(post,bindex(3),bindex(3))));

fprintf(out,' \nAdditional parameters \n');
fprintf(out,'sigma_epsilon^2 , %.4g (%.4g) \n', ...
        mean(chain.sig(post,3,3)-chain.sig(post,3,4)),std(chain.sig(post,3,3)-chain.sig(post,3,4)));
fprintf(out,'sigma_mu^2=sigma_{bar{mu}}^2+sigma_epsilon^2 , %.4g (%.4g) \n', ...
        mean(chain.sig(post,3,3)),std(chain.sig(post,3,3)));
fprintf(out,'rho=corr(mu_lambda_{i1},mu_lambda_{i2}|X,omega,psi) , %.4g (%.4g) \n', ...
        mean(chain.rho(post)),std(chain.rho(post)));
fprintf(out,['corr(mu_lambda_{i1},mu_lambda_{i2}|X) = ' ...
             'sqrt{frac{sigma_{bar{mu}}^2}{sigma_mu^2}} , %.4g (%.4g) \n'], ... 
        mean(sqrt(chain.sig(post,3,4)./chain.sig(post,3,3))), ...
        std(sqrt(chain.sig(post,3,4)./chain.sig(post,3,3))));
fprintf(out,'SigmaKappa, %.4g (%.4g) \n', ...
        mean(chain.sigll(post)),std(chain.sigll(post)));
fprintf(out,'Gamma1, %.4g (%.4g) \n', ...
        mean(chain.theta(post)), std(chain.theta(post)));
fprintf(out,'Gamma2, %.4g (%.4g) \n',...
        mean(chain.shape(post)), std(chain.shape(post)));

%% std errors for implied quantities
options.latentOnly = true;
S = numel(post);
meanlatent = zeros(S,3);
stdlatent = meanlatent;
corrlatent = zeros(S,3,3);
for s=1:S
  for k=1:numel(chain.beta);
    beta{k} =  chain.beta{k}(post(s),:)';
  end
  Sigma = squeeze( chain.sig(post(s),:,:));
  gamma = squeeze( chain.gamma(post(s),:,:));
  theta = squeeze( chain.theta(post(s)));
  shape = squeeze( chain.shape(post(s)));
  rho = squeeze( chain.rho(post(s)));
  bll =  chain.bll(post(s),:)';
  sigll =  chain.sigll(post(s));
  if isfield(chain,'alpha')
    alpha =  chain.alpha(post(s),:)';
  else
    alpha = 0;
  end
  
  setSeed(seed);
  sdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
  meanlatent(s,:) = [mean(sdata.lambda(:)) ...
                     mean(sdata.omega) ...
                     mean(sdata.psi)];
  stdlatent(s,:) =  [std(sdata.lambda(:)) ...
                     std(sdata.omega) ...
                     std(sdata.psi)];
  corrlatent(s,:,:) = [1 mycorr(sdata.eLambda,sdata.omega) ...
                      mycorr(sdata.eLambda,sdata.psi); ...
                      NaN 1 mycorr(sdata.omega,sdata.psi); ...
                      NaN NaN 1];
  latent = {sdata.eLambda, sdata.omega, sdata.psi };
  for j=1:size(simdata.x{1},2)
    for k=1:(numel(latent))
      if (k==1) 
        margCorr(s,k,j) = mycorr(latent{k}(2,:)',sdata.x{3}{2}(:,j));
      else 
        margCorr(s,k,j) = mycorr(latent{k},simdata.x{1}(:,j));
      end
    end
  end  
  
  y = [sdata.eLambda(1,:)'; sdata.eLambda(2,:)'];
  X = [sdata.x{3}{1}; sdata.x{3}{2}];
  rcoef{1}(s,:) = ((X'*X) \ X'*y)';
  for j=2:numel(latent)
    X = sdata.x{j-1};
    rcoef{j}(s,:) = ((X'*X) \ X'*latent{j})';
  end
    
  fprintf('.');
end
fprintf('\n');

fprintf(out,' \nImplied Quantities \n');
fprintf(out,', Lambda , Omega , Psi \n');
fprintf(out,', (Health risk), (Moral hazard) , (Risk aversion)  \n');
fprintf(out,'Expected , %.3g (%.4g) , %.3g  (%.4g) , %.3g  (%.4g)  \n', ...
        mean(simdata.lambda(:)), std(meanlatent(:,1)), ...
        mean(simdata.omega),  std(meanlatent(:,2)), ...
        mean(simdata.psi),  std(meanlatent(:,3)));
fprintf(out,'Std. Dev. , %.3g (%.4g) , %.3g (%.4g) , %.3g (%.4g)  \n', ...
        std(simdata.lambda(:)), std(stdlatent(:,1)), ...
        std(simdata.omega), std(stdlatent(:,2)), ...
        std(simdata.psi), std(stdlatent(:,3)));
fprintf(out,' \nUnconditional correlation matrix  \n');
fprintf(out,', E[Lambda] , Omega , Psi \n');
fprintf(out,' Lambda , %.3g  (%.4g) , %.3g  (%.4g) , %.3g  (%.4g)  \n', ...
        mycorr(simdata.eLambda,simdata.eLambda), std(corrlatent(:,1,1)), ...
        mycorr(simdata.eLambda,simdata.omega), std(corrlatent(:,1,2)), ...
        mycorr(simdata.eLambda,simdata.psi), std(corrlatent(:,1,3)));
fprintf(out,' Omega , -- , %.3g (%.4g) , %.3g (%.4g)  \n', ...
        mycorr(simdata.omega,simdata.omega),std(corrlatent(:,2,2)), ...
        mycorr(simdata.omega,simdata.psi), std(corrlatent(:,2,3)));
fprintf(out,' Psi , -- , -- , %.3g (%.4g)  \n', ...
        mycorr(simdata.psi,simdata.psi), std(corrlatent(:,3,3)));

fprintf(out,' \nMarginal correlations  \n');
fprintf(out,', E[Lambda] , Omega , Psi  \n');
latent = {simdata.eLambda, simdata.omega, simdata.psi };
for j=1:size(simdata.x{1},2)
  fprintf(out,' %s ', xnames{1}{j});
  for k=1:(numel(latent))
    if (k==1) 
      fprintf(out,', %.4g (%.4g)', mycorr(latent{k}(2,:)',simdata.x{3}{2}(:,j)), ...
              std(margCorr(:,k,j)));
    else 
      fprintf(out,', %.4g (%.4g)', mycorr(latent{k},simdata.x{1}(:,j)), ...
               std(margCorr(:,k,j)));
    end
  end
  fprintf(out,' \n');
end  

fprintf(out,' \nRegression coefficients  \n');
fprintf(out,', E[Lambda] , Omega , Psi  \n');
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
      fprintf(out,', %.4g (%.4g)', coef{k}(j),std(rcoef{k}(:,j)));
    else 
      fprintf(out,',      ');
    end
  end
  fprintf(out,' \n');
end  
fclose(out);
