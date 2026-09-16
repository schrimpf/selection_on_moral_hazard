%% Gibbs sampler for insurance choice model
function [chain data] = ...
      gibbs(data, beta0, sig0, g0, rho0, theta0, k0, bll0,sll0, alpha0, opt)
  %pid = 7418;
  global verbosity;
  verbsave = verbosity;
  if ~exist('opt') || isempty(opt) 
    % set default options
    opt.burn = 0; % number of iterations to discard
    opt.save = 1000; % number of iterations to save
    opt.skip = 1; % save every skipth iteration
    opt.display = 100; % display info every this many iterations
    opt.savefile = 'gibbs'; 
    opt.saveInterval = 1e100; % save chain every this many iterations
    opt.useU = false;
    opt.hetsked = false;
    opt.multMH = false;
  end
  if (~isfield(opt,'useU'))
    opt.useU = true;
  end
  if (~isfield(opt,'multMH'))
    opt.multMH = false;
  end

  [order xint wint] = gqzero(data.order);
  order = numel(xint);
  xint = xint'*sqrt(2);
  wint = wint/sqrt(pi);
  total = opt.burn + opt.save*opt.skip;
  beta = beta0;
  sig = sig0;
  gamma = g0;
  theta = theta0;
  shape = k0;
  rho = rho0;
  bll = bll0;
  sigll = sll0;  
  alpha = alpha0;
  for j=1:numel(beta)
    chain.beta{j} = zeros(opt.save,numel(beta{j}));
  end
  chain.sig = zeros(opt.save,size(sig,1),size(sig,2));
  chain.gamma = zeros(opt.save,numel(gamma));
  chain.theta = zeros(opt.save,1);
  chain.shape = zeros(opt.save,1);
  chain.rho = zeros(opt.save,1);
  chain.bll = zeros(opt.save,size(data.xll,2));
  chain.sigll = zeros(opt.save,1);
  % heteroskedasticity : V_i = (1+z_i*alpha).^2 * Sigma
  chain.alpha = zeros(opt.save,numel(alpha));  
  
  nsavelatent = ceil(opt.save/10)+1;
  chain.usave = zeros(nsavelatent,data.N,size(sig,1));
  chain.loglambdaSave = zeros(nsavelatent,data.T,data.N);
  chain.lamloSave = zeros(nsavelatent,data.N);
  nlsaved = 0;
  if (opt.multMH)
    bnd = findOmegaBounds(data);
  else
    bnd = [];
  end
    
  s = 1;
  %u = sampleLatent([],mu,sig,data);
  data.deduct(isnan(data.avail) | data.avail~=true) = -999; 
  % create X for sampling beta
  k = [size(data.x{1},2) size(data.x{2},2) size(data.x{3}{1},2)];
  n = size(data.x{1},1);
  data.X = [ data.x{1}    zeros(n,k(2)) zeros(n,k(3)); ...
             zeros(n,k(1))   data.x{2}   zeros(n,k(3))];
  for t=1:data.T
    data.X = [data.X; zeros(n,k(1)+k(2)) data.x{3}{t}];
  end
  
  %% initial values -- might cause trouble since not rationalizing
  %choices.  
  if (isfield(opt,'restartFile')) 
    load(opt.restartFile);
    s = find(all(chain.beta{1}==0,2),1)-1;
    if (isempty(s))
      s = size(chain.beta{1},1);
    end
    onet = (s-1)*opt.skip+1;
    myprintf('starting with iteration %d, save place %d\n',onet,s);
    nfail = zeros(data.N,1);
  else 
    if (isfield(data,'muL') && opt.useU)
      myprintf(['Using true simulated latent variables as initial ' ...
               'values.\n']);   
      u = [log([data.omega data.psi]) data.muL];
      sigL = data.sigL;  
      loglambda = log(data.lambda+ones(data.T,1)*data.lamlo');
      nfail = zeros(data.N,1);
      lamlo = data.lamlo;
    else % new latent vars
      mu = [data.x{1}*beta{1} data.x{2}*beta{2}];
      for t=1:data.T
        mu = [mu data.x{3}{t}*beta{3}];
      end  
      if (isOctave)
        [r p] = chol_c(sig);
      else
        [r p] = chol(sig);
      end
      u = randn(data.N,2+data.T)*r + mu;
      sigL = 1./sqrt(gamrnd(shape,theta,data.N,1));
      loglambda = (randn(data.N,data.T).*(sigL*ones(1,data.T)) + ...
                   u(:,3:(2+ data.T)))';
      lamlo = randn(data.N,1)*sigll + data.xll*bll;
      % make sure they rationalize choices
      [u sigL lamlo fail] = findValidLatent(u,sigL, lamlo, data, xint, wint,1:data.N, mu, ...
                                            sig,true,opt.multMH);
      nfail = fail;
    end
    %sigL = data.sigL;
    %loglambda = log(data.lambda);
    rootprec = ones(data.N,1); % heteroskedasticity
    
    % save initial values
    for j=1:numel(beta)
      chain.beta{j}(s,:) = beta{j};
    end
    chain.sig(s,:,:) = sig;
    chain.gamma(s,:) = gamma;
    chain.theta(s) = theta;
    chain.shape(s) = shape;
    chain.covu(s,:,:) = cov(u);
    chain.meanu(s,:) = mean(u);
    chain.quantu(s,:,:) = quantile(u,0:0.05:1)';
    chain.meanL(s) = mean(loglambda(:));
    chain.sigL(s) = std(loglambda(:));
    chain.alpha(s,:) = alpha;
    nlsaved = nlsaved+1;
    chain.usave(nlsaved,:,:) = u;
    chain.lamloSave(nlsaved,:) = lamlo;
    chain.loglambdaSave(nlsaved,:,:) = loglambda;
    ndropped = 0;
    onet=1;
  end
  maxFail = 0;
  if (any(nfail>maxFail))
    % drop observations that seems to be impossible to get valid
    ndropped = ndropped + sum(nfail>maxFail);
    fprintf(['WARNING: dropping %d observations because of failure to ' ...
             'find latent variables.  Total dropped=%d\n'],sum(nfail>maxFail),ndropped);
    data = dropBad(data,nfail>maxFail);
    if (opt.multMH)
      bnd = findOmegaBounds(data);
    else
      bnd = [];
    end
    u = u(nfail<=maxFail,:);
    chain.usave = chain.usave(:,nfail<=maxFail,:);
    chain.lamloSave = chain.lamloSave(:,nfail<=maxFail);
    chain.loglambdaSave = chain.loglambdaSave(:,:,nfail<=maxFail);
    sigL = sigL(nfail<=maxFail);
    lamlo = lamlo(nfail<=maxFail);
    loglambda = loglambda(:,nfail<=maxFail);
    rootprec = rootprec(nfail<=maxFail);
    nfail = nfail(nfail<=maxFail);
    if (~opt.multMH)
      [c2 ev]= findChoices_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                             sigL, exp(u(:,2)), ...
                             xint,wint, data.deduct, ...
                             data.maxoop, data.prem, data.avail, ...
                             data.choice, data.totalSpend, ...
                             data.maxIter);
    else 
      [c2 ev]= findChoicesMultMH_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                   sigL, exp(u(:,2)), ...
                                   xint,wint, data.deduct, ...
                                   data.maxoop, data.prem, data.avail, ...
                                   data.choice, data.totalSpend, ...
                                   data.maxIter);
    end

    bad = any(c2~=data.choice & data.choice>0,1);
  end
  if (~opt.multMH)
    [c2 v] = findChoices_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                           sigL, exp(u(:,2)), ...
                           xint,wint, data.deduct, ...
                           data.maxoop, data.prem, data.avail, ...
                           data.choice, data.totalSpend, ...
                           data.maxIter);
  else 
    [c2 v] = findChoicesMultMH_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                 sigL, exp(u(:,2)), ...
                                 xint,wint, data.deduct, ...
                                 data.maxoop, data.prem, data.avail, ...
                                 data.choice, data.totalSpend, ...
                                 data.maxIter);
  end    
  b = (c2~=data.choice & data.choice>0);
  bad = any(c2~=data.choice & data.choice>0,1);
  fprintf('%d bad, %g %g\n',sum(bad),mean(data.choice(find(b))), ...
          mean(c2(find(b))));
  tic
  for t=onet:total % main loop
    if verbsave>1 || (verbsave>0 && (mod(t-1,opt.display)==0 || t<100) )
      verbosity = verbsave;
    else
      verbosity = 0;
    end
      
    if any(isinf(loglambda(:)))
      fprintf('%d loglambda=inf\n',sum(isinf(loglambda(:))));
    end
    %% sample beta
    if (isfield(data,'logOmegaHi') && isfinite(data.logOmegaHi))
      mu = [data.x{1}*beta{1} data.x{2}*beta{2}];
      for t1=1:data.T
        mu = [mu data.x{3}{t1}*beta{3}];
      end  
      if (isOctave)
        [cs p]= chol_c(sig);
      else
        [cs p] = chol(sig);
      end
      e = u;
      e(:,1) = (u(:,1) - mu(:,1))/cs(1,1);
      e(:,2) = (u(:,2) - mu(:,2) - e(:,1)*cs(1,2))/cs(2,2);
      for t1=1:data.T
        e(:,2+t1) = (u(:,2+t1) - mu(:,2+t1) - ...
                    e(:,1:(2+t1-1))*cs(1:(2+t1-1),2+t1))/cs(2+t1,2+t1);
      end
      e = e.*(rootprec*ones(1,size(e,2)));
      unif = normcdf(e(:,1))./ ...
             normcdf((data.logOmegaHi - mu(:,1)).*rootprec/cs(1,1));
      unif(unif>=1) = max(unif(unif<1)); % can happen when searching for
                                         % valid latent
      etilde = e;
      etilde(:,1) = norminv(unif);
      utilde = mu + (etilde*cs)./(rootprec*ones(1,size(u,2)));

      beta = sampleBeta(utilde,sig,data,rootprec);
    else 
      e = u;
      beta = sampleBeta(u,sig,data,rootprec);
    end
    needsig=true;
    tries=0;
    %% sample sigma
    while (needsig)
      if (isfield(data,'logOmegaHi') && isfinite(data.logOmegaHi))
        mu = [data.x{1}*beta{1} data.x{2}*beta{2}];
        for t1=1:data.T
          mu = [mu data.x{3}{t1}*beta{3}];
        end  
        if (isOctave)
          [cs p] = chol_c(sig);
        else 
          [cs p] = chol(sig);
        end
        e = u;     
        e(:,1) = (u(:,1) - mu(:,1))/cs(1,1);
        e(:,2) = (u(:,2) - mu(:,2) - e(:,1)*cs(1,2))/cs(2,2);
        for t1=1:data.T
          e(:,2+t1) = (u(:,2+t1) - mu(:,2+t1) - ...
                      e(:,1:(2+t1-1))*cs(1:(2+t1-1),2+t1))/cs(2+t1,2+t1);
        end
        e = e.*(rootprec*ones(1,size(e,2)));
        unif = normcdf(e(:,1))./ ...
               normcdf((data.logOmegaHi - mu(:,1)).*rootprec/cs(1,1));
        unif(unif>=1) = max(unif(unif<1));
        etilde = e;
        etilde(:,1) = norminv(unif);
        utilde = mu + (etilde*cs)./(rootprec*ones(1,size(etilde,2)));
      else 
        utilde = u;
      end
      e = zeros(size(utilde));
      for j=1:(numel(beta)-1)
        e(:,j) = utilde(:,j) - data.x{j}*beta{j};
      end
      j = numel(beta);
      for t1=1:data.T
        e(:,j+t1-1) = utilde(:,j+t1-1) - data.x{j}{t1}*beta{j};
      end
      e = e.*(rootprec*ones(1,size(e,2)));
      sig(1:2,1:2) = sampleSigmaOP(e,data);
      try 
        gamma(1:2) = sampleGammaOP(e,gamma(3),rho,data);
      catch
        save gibbsCrash.odat;
        error('sampleGammaOp failed');
      end
      gamma(3) = sampleGammaEps(e,gamma(1:2),rho,data);
      rho = sampleRho(e,gamma,rho,data) ;
      sig = combineSigOPbeta(sig(1:2,1:2),gamma,rho,data.T);
      needsig = det(sig)<=0 || ~isreal(sig);
      tries = tries+1;
      %if (~isreal(sig))
      %  error('sig is not real');        
      %end
      if (tries>1)
        myprintf('+');
      end
      if (tries>100)
        error('failed to sample sig 100 times');
      end
    end
    if (tries>1) 
      myprintf('\n');
    end
    %% sample scale and shape of sigL
    theta = sampleTheta(theta,shape,sigL,data);
    shape = sampleShape(shape,theta,sigL,data);
    %% sample bll and sigll
    e = rootprec.*(lamlo - data.xll*bll)/sigll;
    ehi = rootprec.*(data.lamloHi - data.xll*bll)/sigll;
    p = normcdf(e)./normcdf(ehi);
    p(p>1-1e-8) = 1-1e-8;
    p(p<1e-8) = 1e-8;
    ll = data.xll*bll+sigll*norminv(p)./rootprec;
    blllast = bll;
    bll = samplebll(ll,sigll,data,rootprec);
    myprintf('mean xbll = %g\n',mean(data.xll*bll));
    siglllast = sigll;
    sigll = sampleSigll(ll,bll,data,rootprec);
    
    if (any(isnan([beta{1}(:); beta{2}(:); beta{3}(:); sig(:); theta; ...
                   shape; bll(:); sigll])))
      fprintf(['NaN encountered.  Saving current state in tmp file and ' ...
               'stopping.\n']);
      if (isOctave)
        save('-v7',opt.savefile,'chain','data','beta','sig','theta','shape','u','sigL', ... 
             'loglambda','lamlo','bll','sigll','p','ehi','blllast', ...
             'siglllast');
      else 
        save(opt.savefile,'chain','data','beta','sig','theta','shape','u','sigL', ... 
             'loglambda','lamlo','bll','sigll','p','ehi','blllast', ...
             'siglllast');        
      end
      error('NaN encountered.  Saved results in %s',opt.savefile);
    end
    myprintf('before %g %g | %g | %g %g %g %g|\n',mean(u(:,1)), ...
            mean(loglambda(:)),var(u(:,1)),cov(loglambda')) ;
    if ~isequal(size(sigL),size(rootprec))
      sigL = sigL';
    end
    %check(u,loglambda,lamlo,data);
    [u(:,1) loglambda] = ...
        sampleLambdaOmega(u, loglambda, sigL./rootprec, beta, sig, ...
                          data, xint, wint, lamlo, ...
                          rootprec,bnd);    
    if any(isinf(loglambda(:)))
      fprintf('after sampleLambdaOmega %d loglambda=inf\n' ...
              ,sum(isinf(loglambda(:))));
    end
    %check(u,loglambda,lamlo,data);
    myprintf('after %g %g | %g | %g %g %g %g|\n',mean(u(:,1)), ...
            mean(loglambda(:)),var(u(:,1)),cov(loglambda')) ;
    [u(:,2)] = samplePsi(u, loglambda, sigL, beta, sig, ...
                         data, xint, wint, lamlo, rootprec, opt.multMH);

    [u(:,3:(2+data.T))] = sampleMu(u, loglambda, sigL, beta, sig, ...
                                    data, xint, wint, lamlo, rootprec, opt.multMH);
    if any(isnan(sigL) | isinf(sigL))
      sigL
      error('sigL isnaninf');
    end    
    sigL = sampleSigL(u, loglambda, sigL, shape, theta, ...
                      data, xint, wint, lamlo, rootprec, opt.multMH);
    if any(isnan(sigL) | isinf(sigL))
      sigL
      error('sigL isnaninf');
    end
    
    [lamlo loglambda] = sampleLamlo(u, loglambda, sigL, bll, sigll, ...
                                    data,xint,wint,lamlo, ...
                                    rootprec, opt.multMH);
    %myprintf('mean lamlo = %g\n',mean(lamlo));
    %check(u,loglambda,lamlo,data);
    if any(isinf(loglambda(:)))
      fprintf('after sampleLamlo %d loglambda=inf\n',sum(isinf(loglambda(:))));
    end

    % sample alpha and update heteroskedasticity (rootprec)
    if (opt.hetsked)
      [alpha rootprec] = sampleAlpha(alpha,loglambda,u,sigL,beta,sig,lamlo, ...
                                     bll,sigll,data);
    end
    
    if (any(bad)) 
      if (~opt.multMH)
        [c2 ev] = findChoices_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                sigL, exp(u(:,2)), ...
                                xint,wint, data.deduct, ...
                                data.maxoop, data.prem, data.avail, ...
                                data.choice, data.totalSpend, ...
                                data.maxIter);
      else 
        [c2 ev] = findChoicesMultMH_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                      sigL, exp(u(:,2)), ...
                                      xint,wint, data.deduct, ...
                                      data.maxoop, data.prem, data.avail, ...
                                      data.choice, data.totalSpend, ...
                                      data.maxIter);
      end
        
      bad = any(c2~=data.choice & data.choice>0,1);      
      good = find(all(c2==data.choice | data.choice<=0,1));
      fprintf('warning: %d observations are bad\n',sum(bad));      
      %find(bad)
      %data.choice(:,bad)
      %c2(:,bad)
      %ev(:,:,bad)
      [u sigL lamlo fail] = findValidLatent(u,sigL,lamlo, data, xint, wint, find(bad), beta, ...
                                      sig,true,opt.multMH);
      sum(fail)
      nfail(fail) = nfail(fail)+1;
      maxFail = 25;
      if (any(nfail>maxFail))
        % drop observations that seems to be impossible to get valid
        % latent variables for
        ndropped = ndropped + sum(nfail>maxFail);
        fprintf(['WARNING: dropping %d observations because of failure to ' ...
                 'find latent variables.  Total dropped=%d\n'],sum(nfail>maxFail),ndropped);
        data = dropBad(data,nfail>maxFail);
        if (opt.multMH)
          bnd = findOmegaBounds(data);
        else
          bnd = [];
        end
        u = u(nfail<=maxFail,:);
        chain.usave = chain.usave(:,nfail<=maxFail,:);
        chain.lamloSave = chain.lamloSave(:,nfail<=maxFail);
        chain.loglambdaSave = chain.loglambdaSave(:,:,nfail<=maxFail);
        sigL = sigL(nfail<=maxFail);
        lamlo = lamlo(nfail<=maxFail);
        loglambda = loglambda(:,nfail<=maxFail);
        rootprec = rootprec(nfail<=maxFail);
        nfail = nfail(nfail<=maxFail);
        if (~opt.multMH)
          [c2 ev]= findChoices_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                 sigL, exp(u(:,2)), ...
                                 xint,wint, data.deduct, ...
                                 data.maxoop, data.prem, data.avail, ...
                                 data.choice, data.totalSpend, ...
                                 data.maxIter);
        else 
          [c2 ev]= findChoicesMultMH_c(u(:,1),loglambda, lamlo, 1, u(:,3:(2+data.T))', ...
                                       sigL, exp(u(:,2)), ...
                                       xint,wint, data.deduct, ...
                                       data.maxoop, data.prem, data.avail, ...
                                       data.choice, data.totalSpend, ...
                                       data.maxIter);
        end
        bad = any(c2~=data.choice & data.choice>0,1);
      end
    end
    %fprintf('done find\n');
    %system(sprintf('pmap %d | egrep "total|paul"', pid));
    if (t>opt.burn && mod(t-opt.burn-1,opt.skip)==0)
      s = s+1;
      for j=1:numel(beta)
        chain.beta{j}(s,:) = beta{j};
      end      
      chain.sig(s,:,:) = sig;
      chain.gamma(s,:) = gamma;
      chain.rho(s) = rho;
      chain.theta(s) = theta;
      chain.shape(s) = shape;
      chain.bll(s,:) = bll;
      chain.sigll(s) = sigll;
      chain.quantu(s,:,:) = quantile(u,0:0.05:1)';
      chain.meanL(s) = mean(loglambda(:));
      chain.sigL(s) = std(loglambda(:));
      chain.alpha(s,:) = alpha;
      if (mod(s-1,10)==0) 
        nlsaved = nlsaved+1;
        chain.usave(nlsaved,:,:) = u;
        chain.lamloSave(nlsaved,:) = lamlo;
        chain.loglambdaSave(nlsaved,:,:) = loglambda;
      end
    end
    if (mod(t-1,opt.display)==0)      
      toc
      fprintf('\n%d of %d iterations\n',t,total);
      for j=1:numel(beta)
        myprintf(' b{%d} = ',j);
        myprintf('%5.3g ',beta{j});
        myprintf('\n');
      end      
      myprintf('  m(u) = %5.3g %5.3g %5.3g %5.3g\nSigma = \n',mean(u));
      myprintf('   %5.3g %5.3g %5.3g %5.3g \n',sig);
      myprintf('\ncov(u) = \n');
      myprintf('   %5.3g %5.3g %5.3g %5.3g \n',cov(u));
      myprintf(' gamma = %5.3g %5.3g %5.3g \n',gamma);
      myprintf(' rho = %5.3g \n',rho);
      myprintf('ml = %5.3g \n',mean(loglambda(:)));
      myprintf('theta=%5.3g shape=%5.3g\n',theta,shape);
      myprintf('mean(lamlo) = %5.3g\n',mean(lamlo));
      myprintf('alpha = '); myprintf('%5.3g ',alpha); myprintf('\n');
      myprintf('bll = '); myprintf('%5.3g ',bll); myprintf('\n');
      myprintf('\n sigll=%5.3g\n',sigll);
      %if (mod(t,10)==0) 
      %pause();
      %end        
    else % elif !display
      fprintf('.');
    end % fi display
    if (mod(t,opt.saveInterval)==0)
      fprintf('saving intermediary results\n');
      check(u,loglambda,lamlo,data);
      if (isOctave)
        save('-v7',opt.savefile,'chain','data','u','sigL','loglambda', ...
             'lamlo','rootprec','alpha');
      else
        save(opt.savefile,'chain','data','u','sigL','loglambda', ...
             'lamlo','rootprec','alpha');
      end
    end
  end % for(t)
end % function gibbs

%% Given latent variables and Sigma, sample beta
function beta=sampleBeta(u,Sig,data,rootprec);
  if (isempty(data.prior.beta) || any(isinf(data.prior.beta)))
    U = reshape(u,(2+data.T)*data.N,1);    
    if (nargin<4)
      % disperse prior
      isig = Sig \ eye(length(Sig));
      xiSx = XtTimesKronSIdTimesY(data.X,isig,data.N,data.X);
      %oldxiSx = data.X'*kronEye(isig,(data.N))*data.X;    
      v = (xiSx) \ eye(size(data.X,2));
      Bhat = v*XtTimesKronSIdTimesY(data.X,isig,data.N,U); 
    else 
      R = reshape(rootprec*ones(1,size(u,2)),(2+data.T)*data.N,1);
      XR = data.X.*(R*ones(1,size(data.X,2)));
      isig = Sig \ eye(length(Sig));
      xiSx = XtTimesKronSIdTimesY(XR,isig,data.N,XR);
      %oldxiSx = data.X'*kronEye(isig,(data.N))*data.X;    
      v = (xiSx) \ eye(size(data.X,2));
      Bhat = v*XtTimesKronSIdTimesY(XR,isig,data.N,U.*R);       
    end
    if (isOctave)
      [R info] = chol_c(v);
    else
      [R info] = chol(v);
    end
    if (info~=0) 
      fprintf(['WARNING (sampleBeta): v is singular (info=%d), adding ' ...
               'eye(size(v))*norm(v)*1e-6\n'],info);
      R = chol_c(v+eye(size(v))*norm(v)*1e-6);
    end
    B = Bhat + (randn(1,numel(Bhat)) * R)';
    k = 1;
    for j=1:(numel(data.x)-1)
      beta{j} = B(k:(k+size(data.x{j},2)-1));      
      k = k + numel(beta{j});
    end
    j = numel(data.x);
    beta{j} = B(k:(k+size(data.x{j}{1},2)-1));
    k = k + numel(beta{j});
    assert(k==numel(B)+1);
  else 
    error('only disperse prior implemented');
  end
end

%% Given latent variables and beta, sample Sigma
function Sig = sampleSigmaOP(e,data)
  lop = 1:2;
  n = size(e,1);
  nshat = e(:,lop)'*e(:,lop);
  if (~isempty(data.prior.Sigma)) 
    A = nshat + data.prior.Sigma.A;
    m = n + data.prior.Sigma.m;
  else 
    A = nshat;
    m = n;
  end
  Sig = iwishrnd(A,m);
end


%% Sample gammaOP -- assuming diffuse prior
function bop = sampleGammaOP(e,beps, rho, data)
  T = size(e,2)-2;
  lmu = 3:(2+T);
  lop = 1:2;
  n = size(e,1);
  y = reshape(e(:,lmu),T*n,1);
  x = repmat(e(:,lop),T,1);
  ixx = (x'*x)\eye(size(x,2));
  if (isOctave)
    [r p] = chol_c(ixx*beps*beps/(1-rho^2));
  else
    [r p] = chol(ixx*beps*beps/(1-rho^2));
  end
  if (p>0) 
    save problems e beps rho data;
    fprintf('WARNING: sampleGammaOP ixx not positive definite.\n');
    ixx
    beps
    rho
  end
  bop = randn(1,2)*r + (ixx*(x'*y))';
end

%% Sample gammaEps -- assuming diffuse prior
function beps = sampleGammaEps(e, bop, rho, data)
  T = size(e,2)-2;
  lmu = 3:(2+T);
  lop = 1:2;
  n = size(e,1);
  y = reshape(e(:,lmu),data.T*data.N,1);
  x = repmat(e(:,lop),data.T,1);
  beps = 1/sqrt(gamrnd(n*data.T/2, 2/((1-rho^2)*sum((y - x*bop').^2)) ));
end

%% Sample rho -- assuming diffuse prior
function rho = sampleRho(e,gamma,rho,data) 
  T = size(e,2)-2;
  lmu = 3:(2+T);
  lop = 1:2;
  n = size(e,1);
  epsi = (e(:,lmu) - (e(:,lop)*gamma(1:2)')*ones(1,T))/gamma(3);
  %   eem1 = sum(sum(epsi(:,2:T).*epsi(:,1:(T-1))));
  %   em1sq = sum(sum(epsi(:,1:(T-1)).^2));
  %   e1sq = sum(epsi(:,1).^2);
  %   rhoHat = (e1sq + 2*eem1)/(2*em1sq);
  %   rhoTilde = eem1/em1sq;
  %   srho = 1/sqrt(em1sq);
  %rho = randdtn_c(rhoHat,srho,-1,1);
  logf = @(r) sum(sum(-0.5*(epsi(:,2:T) - r*epsi(:,1:(T-1))).^2))+ ...
         sum(-0.5*epsi(:,1).^2*(1-r^2)) + n/2*log(1-r^2);
  logf0 = logf(rho);
  a = 0;
  Tries = 20;
  for t=1:Tries
    rhoT = randdtn_c(rho,1/sqrt(n),-1,1);
    logfTry = logf(rhoT);
    if (log(rand()) < logfTry - logf0)
      logf0 = logfTry;
      rho = rhoT;
      a = a+1;
    end
  end
  myprintf('Accepted %d of %d draws of rho=%g\n',a,Tries,rho);
end

%% Sample Lambda and Omega
function [logomega loglambda] = ...
	  sampleLambdaOmega(u0, logl0, sigL, ...
			    beta,  Sig, data, xint, wint, lamlo, ...
                            rootprec, bnd)
  if (isfield(data,'logOmegaHi'))
    logOmegaHi = data.logOmegaHi;
  else
    logOmegaHi = inf;
  end
  lo = 1;
  lpm = 2:size(u0,2);
  logo0 = u0(:,1);
  if (size(logl0,1)==data.N) 
    logl0 = logl0';
  end
  % conditional mean and std dev of omega
  delta = Sig(lpm,lpm)\Sig(lpm,lo);
  muOm = data.x{1}*beta{1} + (u0(:,2) -data.x{2}*beta{2})*delta(1);
  for t=1:numel(data.x{3}) 
    muOm = muOm + (u0(:,2+t)-data.x{3}{t}*beta{3})*delta(1+t);
  end
  sigOm = sqrt(Sig(lo,lo) - Sig(lo,lpm)*delta);   
  if (isempty(bnd)) % additive model
    [logomega loglambda] = ...
        sampleLambdaOmega_c(logo0,logl0, muOm, sigOm, u0(:,3:size(u0,2))', ...
                            sigL, exp(u0(:,2)), ...
                            xint,wint, data.deduct, ...
                            data.maxoop, data.prem, data.avail, ...
                            data.choice, data.totalSpend, ...
                            data.maxIter, data.threads, ...
                            logOmegaHi, lamlo, rootprec);    
  else % multiplicative model
    [logomega loglambda] =  ...
        sampleLambdaOmegaMultMH_c(logo0,logl0, muOm, sigOm, ...
                                  u0(:,3:size(u0,2))', ...
                                  sigL, exp(u0(:,2)), ...
                                  xint,wint, data.deduct, ...
                                  data.maxoop, data.prem, data.avail, ...
                                  data.choice, data.totalSpend, ...
                                  5, data.threads, ... % 5 = mh reps
                                  logOmegaHi, lamlo, rootprec, bnd);
  end
  %quantile(normcdf(logomega', muOm, sigOm), 0.1:0.1:.9)
end

%% sample theta
function theta = sampleTheta(theta0,shape,sigL,data)
  k = numel(sigL)*shape;
  scale = 1/sum(sigL.^(-2));
  logf = @(th) - data.N*log(1-gamcdf(1/data.varHi,shape,th));
  theta = theta0;
  logf0 = logf(theta);
  Tries = 5;
  a=0;
  for t=1:Tries
    thetaT = 1/gamrnd(k,scale);
    logfT = logf(thetaT);
    if (rand()<exp(logfT-logf0))
      a = a+1;
      theta = thetaT;
      logf0 = logfT;
    end
  end
  myprintf('theta=%g, accepted %d of %d\n',theta,a,Tries);
end

%% Sample psi 
function [logpsi] = samplePsi(u0, loglambda, sigL, beta, Sig, data, xint, ...
                              wint, lamlo, rootprec, multMH)
  % conditional mean and variance 
  T = data.T;
  lo = 1;
  lp = 2;
  lm=3:(2+T);
  mupsi = zeros(1,data.N);
  mu = [data.x{1}*beta{1} data.x{2}*beta{2}];
  for t=1:data.T
    mu=[mu data.x{3}{t}*beta{3}];
  end
  delta = Sig([lo lm],[lo lm]) \ Sig([lo lm],lp);
  mupsi = mu(:,lp) + (u0(:,[lo lm]) - mu(:,[lo lm]))*delta;
  sigpsi = sqrt(Sig(lp,lp) - Sig(lp,[lo lm])*delta);
  logpsi = u0(:,2);
  muL = u0(:,3:(2+T))';
  %myprintf('mean(mupsi)=%g var(mupsi)=%g sigpsi=%g | mean(sigL)=%g var(sigL)=%g\n' ...
  %        ,mean(mupsi), var(mupsi),sigpsi,mean(sigL),var(sigL));
  %myprintf('mean(logpsi)=%g var(logpsi)=%g, corr(sigL,logpsi)=%g\n', ...
  %         mean(logpsi),var(logpsi), corr(logpsi(:),sigL(:)));
  logpsi = samplePsi_c(u0(:,lo),loglambda, mupsi, sigpsi, muL, ...
                       sigL, logpsi, ...
                       xint,wint, data.deduct, ...
                       data.maxoop, data.prem, data.avail, ...
                       data.choice, lamlo, rootprec, ...
                       data.maxIter,data.threads, multMH)';
  %myprintf('after sampling, mean(logpsi)=%g var(logpsi)=%g, corr(sigL,logpsi)=%g\n', ...
  %        mean(logpsi),var(logpsi), corr(logpsi(:),sigL(:)));
end

%% Sample mu
function [muL] = sampleMu(u0, loglambda, sigL, beta, Sig, data, xint, wint, ...
                          lamlo, rootprec, multMH)
  % conditional mean and variance 
  T = data.T;
  lo = 1;
  lp = 2;
  lm=3:(2+T);
  xbO = data.x{1}*beta{1};
  xbP = data.x{2}*beta{2};
  xbL = zeros(T,data.N);
  for t=1:T
    xbL(t,:) = data.x{3}{t}*beta{3};
  end
  muL = u0(:,3:(2+T))';
  muL = sampleMu_c(u0(:,lo),loglambda,u0(:,lp),muL, xbO,xbP,xbL, ...
                   sigL, Sig,...
                   xint,wint, data.deduct, ...
                   data.maxoop, data.prem, data.avail, ...
                   data.choice, lamlo, rootprec, ...
                   data.maxIter,data.threads, multMH)';
end

%% Sample sigL
function sigL = sampleSigL(u0, loglambda, sigL, shape, theta, data, xint, ...
                           wint, lamlo,rootprec, multMH)
  muL = u0(:,(3:size(u0,2)))';
  inc = data.choice>0;
  ki = shape + sum(inc,1)/2;
  thetai = zeros(size(ki));
  if (data.T>1)
    for i=1:data.N
      thetai(i) = 2*theta./(2+theta*sum(( (loglambda(inc(:,i),i)-muL(inc(:,i),i)) ...
                                          .*rootprec(i)).^2));
    end
  else 
    thetai = 2*theta./(2+theta*((loglambda'-muL).*rootprec).^2);
  end
  sigL = sampleSigL_c(u0(:,1),loglambda, ki, thetai, muL, ...
                      sigL, exp(u0(:,2)), ...
                      xint,wint, data.deduct, ...
                      data.maxoop, data.prem, data.avail, ...
                      data.choice, lamlo, ...
                      data.maxIter,data.threads,data.varLo,data.varHi, multMH);
end

%% find values of u and sigL that rationalize choices
function [u sigL lamlo fail] = findValidLatent(u0,sigL0, lamlo, data, xint, wint, badindex, ...
                                               mu, sig, newbounds, multMH)
  persistent lb;
  persistent ub;
  if ~exist('newbounds')
    newbounds = false;
  end
  if (size(sigL0,1)==1)
    sigL0 =sigL0';
  end
  if ~exist('multMH')
    multMH = false;
  end
  if (isempty(lb) || newbounds)
    N = data.N;
    T = data.T;
    lb = ones(N,1)*mean([u0 log(sigL0)]) - (3+rand(N,1))*std([u0 log(sigL0)]);
    ub = ones(N,1)*mean([u0 log(sigL0)]) + (3+rand(N,1))*std([u0 ...
                        log(sigL0)]);
    for t=1:T
      ub(ub(:,2+t)>log(1e5),2+t)=log(1e5);
    end
    lohi = min(data.logOmegaHi,log(5e4));
    ub(ub(:,1)>lohi ,1) = lohi;
    lb(lb(:,end)<log(data.varLo),end)=log(data.varLo);
    ub(ub(:,end)>log(data.varHi),end)=log(data.varHi);
    lb = [lb zeros(N,1)];
    ub = [ub ones(N,1)*data.lamloHi];
    lb(lb(:,2)<-20,2)=-20;
    lb=lb';
    ub=ub';
  end
  u = u0;
  sigL = sigL0;
  muL = u0(:,3:size(u0,2))';
  tic;
  if (~multMH) 
    [u sigL lamlo fail]=findValidLatent3_c(muL, log(sigL0), u0(:,1), u0(:,2), ...
                                           data.deduct, data.maxoop, ...
                                           data.prem, xint, wint, data.choice, ...
                                           data.threads, ...
                                           lb,ub, lamlo,data.totalSpend);  
  else
    [u sigL lamlo fail]=findValidLatentMultMH_c(muL, log(sigL0), u0(:,1), u0(:,2), ...
                                                data.deduct, data.maxoop, ...
                                                data.prem, xint, wint, data.choice, ...
                                                data.threads, ...
                                                lb,ub, lamlo,data.totalSpend);  
  end
  toc;
  u = u';
  [mean(u); mean(u0)]
  [var(u); var(u0)]
  [mean(sigL) mean(sigL0)]
  [var(sigL) var(sigL0)]
  myprintf('Failed %d times\n',sum(fail));
  %pause();
  
  %   [ev ul] = expectedValue(u(~fail,:),sigL(~fail), xint,wint, ...
  %                           squeeze(data.deduct(:,4,~fail)), ...
  %                           squeeze(data.maxoop(:,4,~fail)), ...
  %                           squeeze(data.prem(:,4,~fail)));

end

%% sample shape parameter
function shape = sampleShape(shape,theta,sigL,data)
  shape = sampleShape_c(shape,theta,sigL.^(-2),1/data.varHi);
end

%% bll
function b = samplebll(lamlo,sig,data,rootprec)
  if (nargin>3)
    x = data.xll.*(rootprec*ones(1,size(data.xll,2)));
    y = lamlo.*rootprec;
  else 
    x = data.xll;
    y = lamlo;
  end
  xx = x'*x;
  mb = xx \ (x'*y);
  vb = xx \ eye(size(xx))*sig*sig;
  if (isOctave)
    [r p] = chol_c(vb);
  else
    [r p] = chol(vb);
  end
  b = mb + (randn(1,numel(mb)) * r)';
end

function sig = sampleSigll(lamlo,b,data,rootprec)
  if (nargin>3)
    e = (lamlo - data.xll*b).*rootprec;
  else
    e = lamlo - data.xll*b;
  end
  sig = sqrt(1/gamrnd(numel(lamlo)/2,2/(sum(e.^2))));
end

function [lamlo loglambda] = sampleLamlo(u, loglambda, sigL, bll, sigll,...
                                         data,xint,wint,lamlo, rootprec, multMH)
  T = size(loglambda,1);
  %check(u,loglambda,lamlo,data);
  ll0 = loglambda;
  lam = exp(loglambda) - (lamlo*ones(1,T))';  
  muL = u(:,3:end)';
  mull = data.xll*bll;
  ll0 = loglambda;
  lamlo0 = lamlo;
  [lamlo loglambda] = sampleLamlo_c(lam, muL, sigL, mull, sigll, lamlo0, ...
                                    data.choice, exp(u(:,1)), exp(u(:,2)), ...
                                    xint, wint, data.deduct, data.maxoop, ...
                                    data.prem,data.lamloHi, rootprec, data.threads, ...
                                    multMH);
  %lam0 = lam;
  if any(isinf(loglambda))
    loglambda(isinf(loglambda)) = ll0(isinf(loglambda));
  end
  %lam = exp(loglambda) - (lamlo*ones(1,T))'; 
  %check(u,loglambda,lamlo,data);
  %mean(loglambda(:))
end

function [u sigL lamlo loglambda data] = loadLatent(filename)
  load(filename,'u','sigL','lamlo','loglambda','data');
end

function [alpha rootprec] = sampleAlpha(alpha0,loglambda,u,sigL,Beta,Sig,lamlo, ...
                                        bll,sll,data)
  persistent sigAlpha;
  persistent accept;
  persistent tries;
  if (isempty(sigAlpha))
    sigAlpha = 0.02;
  end
  if (isempty(accept))
    accept = 0;
  end
  if (isempty(tries))
    tries = 0;
  end
  e = zeros(size(u));
  for j=1:2
    e(:,j) = u(:,j) - data.x{j}*Beta{j};
  end
  j = 3;
  for t=1:numel(data.x{3})
    e(:,j) = u(:,j) - data.x{3}{t}*Beta{3};
    j = j+1;
  end
  w = @(a) (1+data.xhs*a).^2;
  et = (loglambda' - e(:,3:end));
  isig = Sig \ eye(length(Sig));
  ell = (lamlo-data.xll*bll)/sll;
  ewe = zeros(data.N,1);
  for i=1:data.N;
    ewe(i) = e(i,:)*isig*e(i,:)' + et(i,:)*et(i,:)'./sigL(i).^2  + ...
             ell(i);
  end
  loglike = @(a) -0.5*sum((ewe)./w(a)) - sum((size(e,2)+data.T+1)* ...
                                             log(w(a)));
  alpha = alpha0;
  for t=1:5
    atry = alpha + sigAlpha * randn(size(alpha));
    u = rand;
    if (u<exp(loglike(atry)-loglike(alpha)))
      alpha = atry;
      accept = accept + 1;
    end
    tries = tries+1;
  end
  rootprec = 1./sqrt(w(alpha));
  myprintf('sampleAlpha: accepted %d of %d, rate=%3.2f\n',accept,tries,accept/tries);
end

%% wrapper for fprintf that hides output when verbosity = 0
function myprintf(varargin)
  global verbosity;
  if (verbosity>0)
    fprintf(varargin{:});
  end
end 