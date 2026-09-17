%% Simulates data with N observations using the given parameters
% Draws insurance plans and covariates from data with replacement
function sim = simulate(N,data,beta,Sig,gamma,rho,k,theta,bll,sll,alpha,options)
  if (isempty(options))
    options = alpha;
    alpha = 0;
  end
  if (~isstruct(options))
    options.panel = true;
  end

  if (~isfield(options,'panel'))
    options.panel = true;
  end
  fieldnames = {'randChoice','sampleCondChoice','constMu', ...
                'cf','sampleAll','avgPlan','nosel','avgX','uncorrX', ...
                'multMH','latentOnly'};
  for f=1:numel(fieldnames)
    if (~isfield(options,fieldnames{f}))
      options.(fieldnames{f}) = false;
    end
  end
  if (~isfield(options,'balance15'))
    options.balance15 = -1;
  end
  if (options.multMH)
    findChoices = @(varargin) findChoicesMultMH_c(varargin{:});
    findChoicesNMH = @(varargin) findChoicesNmultMH_c(varargin{:});
  else
    findChoices = @(varargin) findChoices_c(varargin{:});
    findChoicesNMH = @(varargin) findChoicesNMH_c(varargin{:});
  end

  parm.beta = beta;
  parm.sig=Sig;
  parm.gamma=gamma;
  parm.shape = k;
  parm.theta=theta;
  parm.bll = bll;
  parm.sigll = sll;
  parm.alpha = alpha;

  sim.parm=parm;
  sim.N = N;
  if (options.panel);
    if (options.sampleAll)
      i = kron((1:data.N)',ones(N/data.N,1));
    else
      i = ceil(rand(N,1)*data.N);
    end
    sim.id = (1:N)';
    sim.year = unique(data.year);
    sim.T = numel(unique(sim.year));
    sim.uid = sim.id;
  else
    i = ceil(rand(N,1)*data.N);
    sim.id = (1:N)';
    sim.year = mode(data.year)*size(sim.id);
    sim.T = 1;
  end
  T = sim.T;
  sim.choice = zeros(T,N);
  sim.totalSpend = zeros(T,N);
  % sample observations with replacement
  fields = {'avail','prem','deduct','copay','maxoop'};
  if (~options.avgX)
    if (options.uncorrX)
      for j=1:(numel(data.x)-1)
        i = ceil(rand(N,1)*data.N);
        sim.x{j} = data.x{j}(i,:);
      end
      j = numel(data.x);
      i = ceil(rand(N,1)*data.N);
      for t=1:data.T
        sim.x{j}{t} = data.x{j}{t}(i,:);
      end
      sim.xll = data.xll(i,:);
      sim.xhs = data.xhs(i,:);
    else
      for j=1:(numel(data.x)-1)
        sim.x{j} = data.x{j}(i,:);
      end
      j = numel(data.x);
      for t=1:data.T
        sim.x{j}{t} = data.x{j}{t}(i,:);
      end
      sim.xll = data.xll(i,:);
      sim.xhs = data.xhs(i,:);
    end
  else
    for j=1:(numel(data.x)-1)
      sim.x{j} = ones(sim.N,1)*mean(data.x{j},1);
    end
    j = numel(data.x);
    for t=1:data.T
      sim.x{j}{t} = ones(sim.N,1)*mean(data.x{j}{t},1);
    end
    sim.xll = ones(sim.N,1)*mean(data.xll,1);
    sim.xhs = ones(sim.N,1)*mean(data.xhs,1);
  end
  for field = fields
    if (options.panel)
      sim.(field{1}) = data.(field{1})(:,:,i);
    else
      sim.(field{1}) = data.(field{1})(:,i);
    end
  end
  sim.choiceObs = data.choice(:,i);
  sim.spendObs = data.totalSpend(:,i);
  sim.order = data.order;
  sim.varLo = data.varLo;
  sim.varHi = data.varHi;
  sim.maxIter = data.maxIter;
  sim.threads=data.threads;
  % stuff for taking expectation over lambda given other shocks
  [order xint wint] = gqzero(sim.order);
  xint = xint'*sqrt(2);
  order = numel(xint);
  wint = wint/sqrt(pi);
  % sample shocks
  mu = [sim.x{1}*beta{1} sim.x{2}*beta{2}];
  for t=1:T
    mu = [mu sim.x{3}{t}*beta{3}];
  end

  %% make everyone face the average select plan
  covg = ones(data.N,1);
  scovg = ones(sim.N,1);
  for t=2:4
    covg(data.x{1}(:,t)==1) = t;
    scovg(sim.x{1}(:,t)==1) = t;
  end
  sim.covg = scovg;
  if (options.avgPlan)
    sel = squeeze(data.avail(5,2,:)==1 & data.deduct(5,2,:)>=0);
    fields = {'prem','deduct','maxoop','copay'};
    for t=1:4
      for f=1:numel(fields)
        sim.([fields{f} 'Obs']) = sim.(fields{f});
        sim.(fields{f})(:,2,scovg==t) = mean(data.(fields{f})(:,2,covg==t & sel),3) ...
            * ones(1,sum(scovg==t));
      end
    end
    sim.avail(:,2,:) = 1;
  end
  if ~isempty(alpha)
    sim.stdevi = sqrt((1+sim.xhs*alpha).^2 /  mean((1+sim.xhs*alpha).^2));
  else
    sim.stdevi = ones(size(sim.xll,1),1);
  end

  if (options.sampleCondChoice)
    sim.choice = data.choice(:,i);
    bad = ones(N,1);
    sim.omega = zeros(N,1);
    sim.psi = zeros(N,1);
    sim.muL= zeros(N,T);
    sim.sigL= zeros(N,1);
    tries = 0;
    while (any(bad))
      if (isfield(data,'logOmegaHi'))
        e = randn(N,2+T);
        cs = chol(Sig);
        e(:,1) = norminv(rand(N,1).*normcdf((data.logOmegaHi-mu(:,1))./(sim.stdevi*cs(1,1))));
        assert(all(e(:,1).*sim.stdevi*cs(1,1)+mu(:,1) < data.logOmegaHi),'bad u');
      else
        e = randn(N,2+T);
      end
      u = e.*(sim.stdevi*ones(1,size(e,2)))*chol(Sig) + mu;

      sim.omega(bad) = exp(u(bad,1));
      sim.psi(bad) = exp(u(bad,2));
      sim.muL(bad,:) = u(bad,3:(2+T));
      sim.sigL(bad) = sqrt(1./gamrnd(k,theta,sum(bad),1));
      ll =sim.xll*bll + randn(N,1).*sim.stdevi*sll;
      sim.lamlo(bad) = ll(bad);
      badsl = sim.sigL<sqrt(sim.varLo) | sim.sigL>sqrt(sim.varHi);
      while (any(badsl))
        sim.sigL(badsl) = sqrt(1./gamrnd(k,theta,sum(badsl),1));
        badsl = sim.sigL<sqrt(sim.varLo) | sim.sigL>sqrt(sim.varHi);
      end
      loglambda = [];
      sim.deduct(isnan(sim.deduct))=-999;
      fprintf('finding choices ...\n');
      [c2 sim.value]= findChoices(log(sim.omega),loglambda, sim.lamlo, 1, ...
                                  sim.muL', sim.sigL, sim.psi, ...
                                  xint,wint, sim.deduct, ...
                                  sim.maxoop, sim.prem, sim.avail, ...
                                  sim.choice, sim.totalSpend, ...
                                  sim.maxIter);
      good = (c2 == sim.choice) | sim.choice<1;
      bad = any(~good,1);
      tries = tries+1;
      if (mod(tries,10)==0)
        fprintf('%d: %d bad\n',tries,sum(bad));
      end
      if (tries>100)
        fprintf('too many tries. giving up for %d\n',sum(bad));
        break;
      end
    end
  else
    if (isfield(data,'logOmegaHi'))
      e = randn(N,2+T);
      cs = chol(Sig);
      e(:,1) = norminv(rand(N,1).*normcdf((data.logOmegaHi-mu(:,1))./(sim.stdevi*cs(1,1))));
      assert(all(e(:,1).*sim.stdevi*cs(1,1)+mu(:,1) < data.logOmegaHi),'bad u');
    else
      e = randn(N,2+T);
    end
    u = e.*(sim.stdevi*ones(1,size(e,2)))*chol(Sig) + mu;
    sim.sigL = sqrt(1./gamrnd(k,theta,N,1));
    badsl = sim.sigL<sqrt(sim.varLo) | sim.sigL>sqrt(sim.varHi);
    while (any(badsl))
      sim.sigL(badsl) = sqrt(1./gamrnd(k,theta,sum(badsl),1));
      badsl = sim.sigL<sqrt(sim.varLo) | sim.sigL>sqrt(sim.varHi);
    end
    if (isfield(data,'lamloHi'))
      p = normcdf(data.lamloHi,sim.xll*bll,sll);
      r = rand(N,1).*p;
      sim.lamlo = norminv(r).*sim.stdevi*sll + sim.xll*bll;
      % randdtn_c(sim.xll*bll,sll*ones(N,1),-1e300*ones(N,1),data.lamloHi*ones(N,1));
    else
      sim.lamlo = sim.xll*bll + randn(N,1).*sim.stdevi*sll;
    end
  end
  if (isfield(data,'logOmegaHi'))
    sim.logOmegaHi = data.logOmegaHi;
  end
  sim.omega = exp(u(:,1));
  sim.psi = exp(u(:,2));
  sim.muL = u(:,3:(2+T));
  if (options.constMu)
    for t=2:T
      sim.muL(:,t) = sim.muL(:,1);
    end
  end
  if (size(sim.lamlo,1)==1)
    sim.lamlo=sim.lamlo';
  end

  loglam = sim.muL(:,1)*ones(1,order) + sim.sigL*xint;
  lam = exp(loglam) - sim.lamlo*ones(1,order);
  sim.eLambda(1,:) = (lam*wint)';

  loglam = sim.muL(:,2)*ones(1,order) + sim.sigL*xint;
  lam = exp(loglam) - sim.lamlo*ones(1,order);
  sim.eLambda(2,:) = (lam*wint)';

  sim.lambda = exp(randn(T,N).*((sim.sigL.*sim.stdevi)*ones(1,T))' + ...
                   sim.muL') - (sim.lamlo*ones(1,T))';

  if (options.latentOnly)
    return;
  end

  %sim.avail(:,5) = false;
  if (options.balance15>=0)
    sim=balance15(sim,options.balance15);
  end % options.balance15

  %% predicted choices
  fprintf('finding choices ...\n');
  [c2 sim.value sim.exSpend eo] = ...
      findChoices(log(sim.omega),1, sim.lamlo, 1, sim.muL', ...
                  sim.sigL, sim.psi, ...
                  xint,wint, sim.deduct, ...
                  sim.maxoop, sim.prem, sim.avail, ...
                  sim.choice, sim.totalSpend, ...
                  sim.maxIter);
  fprintf('found choices\n');
  %% random choices
  sim.randChoice = sim.choice;
  if (options.avgPlan)
    flex = squeeze(sim.avail(5,:,:)~=1) & sim.choice>0;
    sel = squeeze(sim.avail(5,:,:)==1) & sim.choice>0;
  else
    flex = squeeze(data.avail(5,:,:)~=1) & data.choice>0;
    sel = squeeze(data.avail(5,:,:)==1) & data.choice>0;
  end
  pf = zeros(1,5);
  ps = zeros(1,5);
  for c=1:5
    if (options.avgPlan || options.balance15>=0)
      if any(flex(:))
        pf(c) = mean(double(sim.choice(flex)==c));
      end
      if any(sel(:))
        ps(c) = mean(double(sim.choice(sel)==c));
      end
    else
      if any(flex(:))
        pf(c) = mean(double(data.choice(flex)==c));
      end
      if any(sel(:))
        ps(c) = mean(double(data.choice(sel)==c));
      end
    end
  end
  sflex = squeeze(sim.avail(5,:,:)~=1);
  ssel = squeeze(sim.avail(5,:,:)==1);
  if sum(pf)>0
    pf = pf / sum(pf);
    if sum(sflex(:))>0
      sim.randChoice(sflex)=sum(mnrnd(1,pf,sum(sflex(:))).*(ones(sum(sflex(:)),1)*(1:5)),2);
    end
  end
  if sum(ps)>0
    ps = ps / sum(ps);
    if sum(ssel(:))>0
      sim.randChoice(ssel)=sum(mnrnd(1,ps,sum(ssel(:))).*(ones(sum(ssel(:)),1)*(1:5)),2);
    end
  end
  sim.randChoice(all(sim.deduct<0 | isnan(sim.deduct),1)) = -1;
  sim.deduct(isnan(sim.deduct))=-999;
  % find chosen plans
  if (options.randChoice)
    sim.choice = sim.randChoice;
  else
    sim.choice = c2;
  end
  sim.choice(all(sim.deduct<0,1)) = -1;

  % compute spending
  [sim.totalSpend sim.oop sim.oopRate] = spending(sim);

  sim.select = squeeze(sim.avail(5,:,:)==1);
  if (options.cf)
    % no moral hazard
    nmdata=sim;
    [sim.spendNMHpchoice sim.oopNMHpchoice] = spending(nmdata,true);
    [nmdata.choice sim.valueNMH] = ...
        findChoicesNMH(log(sim.omega),1, sim.lamlo, 1, sim.muL', ...
                       sim.sigL, sim.psi, ...
                       xint,wint, sim.deduct, ...
                       sim.maxoop, sim.prem, sim.avail, ...
                       sim.choice, sim.totalSpend, ...
                       sim.maxIter);
    fprintf('found choices nmh');
    sim.choiceNMH = nmdata.choice;
    [sim.spendNMH sim.oopNMH] = spending(nmdata,true);
    % random choices
    rcdata = sim;
    rcdata.choice=sim.randChoice;
    [sim.spendRandChoice sim.oopRandChoice] = spending(rcdata);
    [sim.spendRandChoiceNMH sim.oopRandChoiceNMH] = spending(rcdata,true);

    % no selection on omega
    rcdata = sim;
    ri = ceil(rand(N,1)*N);
    if (isfield(data,'logOmegaHi'))
      rcdata.omega = norminv(rand(N,1).* ...
                             normcdf((data.logOmegaHi-mu(ri,1))./(sim.stdevi*sqrt(Sig(1,1)))));
    else
      rcdata.omega = randn(N,1);
    end
    rcdata.omega = exp(rcdata.omega.*sim.stdevi*sqrt(Sig(1,1)) + mu(ri,1));
    [sim.spendNoOmSel sim.oopNoOmSel] = spending(rcdata);
    [sim.spendNoOmSelNMH sim.oopNoOmSelNMH] = spending(rcdata,true);
    sim.omegaNoSel = rcdata.omega;
    [junk sim.valueNoOmSel foo bar] = findChoices(log(rcdata.omega),log(sim.lambda), sim.lamlo, 1, ...
                                                  sim.muL', sim.sigL, sim.psi, ...
                                                  xint,wint, sim.deduct, ...
                                                  sim.maxoop, sim.prem, sim.avail, ...
                                                  sim.choice, sim.totalSpend, ...
                                                  sim.maxIter);
    [junk sim.valueNoOmSelNMH] = ...
        findChoicesNMH(log(rcdata.omega),1, sim.lamlo, 1, sim.muL', ...
                       sim.sigL, sim.psi, ...
                       xint,wint, sim.deduct, ...
                       sim.maxoop, sim.prem, sim.avail, ...
                       sim.choice, sim.totalSpend, ...
                       sim.maxIter);
    fprintf('found choices nmh 2');
    % no select
    noseldata = sim;
    select = squeeze(sim.avail(5,2,:)==1);
    noseldata.prem(:,2,select) = noseldata.prem(:,1,select);
    noseldata.maxoop(:,2,select) = noseldata.maxoop(:,1,select);
    noseldata.deduct(:,2,select) = noseldata.deduct(:,1,select);
    noseldata.avail(4:5,2,:) = 0;
    [noseldata.choice sim.valueNoSel foo bar] = findChoices(log(noseldata.omega),1, noseldata.lamlo, 1, ...
                                                      noseldata.muL', noseldata.sigL, noseldata.psi, ...
                                                      xint,wint, noseldata.deduct, ...
                                                      noseldata.maxoop, noseldata.prem, noseldata.avail, ...
                                                      noseldata.choice, noseldata.totalSpend, ...
                                                      noseldata.maxIter);
    fprintf('found choices 2');
    noseldata.choice(all(noseldata.deduct<0,1)) = -1;
    [noseldata.totalSpend noseldata.oop noseldata.oopRate] = spending(noseldata);
    sim.choiceNoSel = noseldata.choice;
    sim.spendNoSel = noseldata.totalSpend;
    sim.oopNoSel = noseldata.oop;
    % no heterogeneity in omega / mu_lambda
    logomega = log(mean(sim.omega))*ones(size(sim.omega));
    [sim.choiceNoOmHet sim.valueNoOmHet foo bar]= findChoices(logomega,1, sim.lamlo, 1, ...
                                                      sim.muL', sim.sigL, sim.psi, ...
                                                      xint,wint, sim.deduct, ...
                                                      sim.maxoop, sim.prem, sim.avail, ...
                                                      sim.choice, sim.totalSpend, ...
                                                      sim.maxIter);
    fprintf('found choices 3');
    sdata = sim;
    sdata.choice = sim.choiceNoOmHet;
    sdata.omega = exp(logomega);
    [sim.spendNoOmHet sim.oopNoOmHet sim.oopRateNoOmHet] = spending(sdata);
    sdata.choice = sim.randChoice;
    [sim.spendRandChoiceNoOmHet sim.oopRandChoiceNoOmHet] = spending(sdata);
    sim.spendNoOmHetNMH  = spending(sdata,true);

    logomega = mean(log(sim.omega))*ones(size(sim.omega));
    [sim.choiceNoOmHetLog sim.valueNoOmHetLog foo bar]= findChoices(logomega,1, sim.lamlo, 1, ...
                                     sim.muL', sim.sigL, sim.psi, ...
                                     xint,wint, sim.deduct, ...
                                     sim.maxoop, sim.prem, sim.avail, ...
                                     sim.choice, sim.totalSpend, ...
                                     sim.maxIter);
    fprintf('found choices 4');
    sdata = sim;
    sdata.choice = sim.choiceNoOmHetLog;
    sdata.omega = exp(logomega);
    [sim.spendNoOmHetLog sim.oopNoOmHetLog] = spending(sdata);
    sdata.choice = sim.randChoice;
    [sim.spendRandChoiceNoOmHetLog sim.oopRandChoiceNoOmHetLog] = spending(sdata);
    sim.spendNoOmHetLogNMH  = spending(sdata,true);

    muL = ones(sim.N,1)*mean(sim.muL);
    sim.choiceNoMuLHet= findChoices(log(sim.omega),log(sim.lambda), sim.lamlo, 1, ...
                                      muL', sim.sigL, sim.psi, ...
                                      xint,wint, sim.deduct, ...
                                      sim.maxoop, sim.prem, sim.avail, ...
                                      sim.choice, sim.totalSpend, ...
                                      sim.maxIter);
    fprintf('found choices 6');
    sdata = sim;
    sdata.choice = sim.choiceNoMuLHet;
    sdata.muL = muL;
    [sim.spendNoMuLHet sim.oopNoMuLHet] = spending(sdata);
    sim.spendNoMuLHetNMH  = spending(sdata,true);
    fields = {'value','valueNoSel','valueNMH','valueNoOmHet','valueNoOmHetLog', ...
              'valueNoOmSel','valueNoOmSelNMH'};
    for c=1:5
      for f=1:numel(fields)
        sim.(fields{f})(c,:,:) = squeeze(sim.(fields{f})(c,:,:))./ ...
            (ones(sdata.T,1)*sdata.psi');
        sim.(fields{f})(~(sim.avail==1)) = NaN;
      end
    end
  end % if options.cf

  loglam = sim.muL(:,1)*ones(1,order) + sim.sigL*xint;
  lam = exp(loglam) - sim.lamlo*ones(1,order);
  sim.eLambda(1,:) = (lam*wint)';

  loglam = sim.muL(:,2)*ones(1,order) + sim.sigL*xint;
  lam = exp(loglam) - sim.lamlo*ones(1,order);
  sim.eLambda(2,:) = (lam*wint)';
  sim.eSpendNI = max(lam,0)*wint;
  sim.eSpendFI = max(lam + sim.omega*ones(1,order),0)*wint;
  sim.eMoralHaz = sim.eSpendFI - sim.eSpendNI;
  sim.spendFI = max(sim.lambda+ones(sim.T,1)*sim.omega',0);
  sim.spendNI = max(sim.lambda,0);

  if (options.nosel) %% create spending and choices without select
                     %options too
    ns = sim;
    select = squeeze(ns.avail(4,2,:)==1);
    ns.prem(:,2,select) = ns.prem(:,1,select);
    ns.maxoop(:,2,select) = ns.maxoop(:,1,select);
    ns.deduct(:,2,select) = ns.deduct(:,1,select);
    ns.deduct(4:5,2,:) = -999;
    ns.avail(4:5,2,:) = 0;
    ns.choice = findChoices(log(ns.omega),1, ns.lamlo, 1, ...
                              ns.muL', ns.sigL, ns.psi, ...
                              xint,wint, ns.deduct, ...
                              ns.maxoop, ns.prem, ns.avail, ...
                              ns.choice, ns.totalSpend, ...
                              ns.maxIter);
    ns.choice(all(ns.deduct<0,1)) = -1;
    [ns.totalSpend ns.oop ns.oopRate] = spending(ns);
    sim.spendNoSelect = ns.totalSpend;
    sim.oopRateNoSelect = ns.oopRate;
    sim.choiceNoSelect = ns.choice;
  end

end
