function results = test_gibbs_pipeline()
%% test_gibbs_pipeline.m
% Integration test: runs a multi-step Gibbs sampling cycle on a self-contained
% synthetic dataset to verify end-to-end pipeline consistency.

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_gibbs_pipeline (synthetic end-to-end integration)...\n');

try
  % 1. Create synthetic panel data (N = 40, T = 2, K = 3 plans)
  N = 40;
  T = 2;
  np = 3;

  % Design matrices: 2 covariates for omega, 2 for psi, 2 for lambda
  kx = 2;
  data.N = N;
  data.T = T;
  data.x{1} = [ones(N, 1), randn(N, kx-1)];
  data.x{2} = [ones(N, 1), randn(N, kx-1)];
  for t = 1:T
    data.x{3}{t} = [ones(N, 1), randn(N, kx-1)];
  end
  data.X = blkdiag(data.x{1}, data.x{2}, data.x{3}{1}, data.x{3}{2});

  data.xll = [ones(N, 1), randn(N, 1)];
  data.xhs = ones(N, 1);
  data.threads = 4;
  data.maxIter = 100;
  data.varLo = 0.01;
  data.varHi = 10.0;
  data.logOmegaHi = Inf;
  data.lamloHi = 2000;

  % Plan menus
  data.deduct = zeros(np, T, N);
  data.maxoop = zeros(np, T, N);
  data.prem = zeros(np, T, N);
  data.avail = ones(np, T, N);
  for i = 1:N
    data.deduct(:, :, i) = repmat([1000; 500; 250], 1, T);
    data.maxoop(:, :, i) = repmat([3000; 1500; 750], 1, T);
    data.prem(:, :, i)   = repmat([0; 200; 500], 1, T);
  end

  data.totalSpend = 500 * rand(T, N);
  data.choice = ones(T, N);
  data.default = zeros(T, N);

  [~, xint, wint] = gqzero(24);
  xint = xint'*sqrt(2);
  wint = wint/sqrt(pi);

  % 2. Initial parameters
  beta{1} = [log(1000); 0];
  beta{2} = [-5.0; 0];
  beta{3} = [6.0; 0];
  gamma = [0.1; 0.1; 0.5];
  rho = 0.3;
  theta = 1.0;
  shape = 2.0;
  bll = [50; 0];
  sigll = 10.0;
  Sigma = combineSigOPbeta(0.5*eye(2), gamma, rho, T);
  rootprec = ones(N, 1);

  % 3. Initialize latent variables
  mu = [data.x{1}*beta{1}, data.x{2}*beta{2}];
  for t = 1:T
    mu = [mu, data.x{3}{t}*beta{3}];
  end
  [R_sig, ~] = chol_c(Sigma);
  u = randn(N, 2+T)*R_sig + mu;
  sigL = 1./sqrt(gamrnd(shape, theta, N, 1));
  sigL = sigL(:);
  loglambda = (randn(N, T).*(sigL*ones(1, T)) + u(:, 3:(2+T)))';
  lamlo = 50 * ones(N, 1);

  % Rationalize initial choices
  [choice_init, ~] = findChoices_c(u(:,1), loglambda, lamlo, 1, u(:,3:end)', ...
                                  sigL, exp(u(:,2)), xint, wint, data.deduct, ...
                                  data.maxoop, data.prem, data.avail, data.choice, ...
                                  data.totalSpend, data.maxIter);
  data.choice = choice_init;

  % 4. Run 3 Gibbs iterations
  n_iters = 3;
  for iter = 1:n_iters
    % A. sampleBeta
    U = reshape(u, (2+T)*N, 1);
    R = reshape(rootprec*ones(1, size(u,2)), (2+T)*N, 1);
    XR = data.X.*(R*ones(1, size(data.X, 2)));
    isig = Sigma \ eye(length(Sigma));
    xiSx = XtTimesKronSIdTimesY(XR, isig, N, XR);
    v = (xiSx) \ eye(size(data.X, 2));
    Bhat = v*XtTimesKronSIdTimesY(XR, isig, N, U.*R);
    [Rchol, ~] = chol_c(v);
    B = Bhat + (randn(1, numel(Bhat)) * Rchol)';
    k = 1;
    for j = 1:(numel(data.x)-1)
      beta{j} = B(k:(k+size(data.x{j},2)-1));
      k = k + numel(beta{j});
    end
    beta{numel(data.x)} = B(k:(k+size(data.x{numel(data.x)}{1},2)-1));

    % B. sampleSigmaOP
    e = zeros(size(u));
    for j = 1:(numel(beta)-1)
      e(:,j) = u(:,j) - data.x{j}*beta{j};
    end
    j = numel(beta);
    for t1 = 1:T
      e(:,j+t1-1) = u(:,j+t1-1) - data.x{j}{t1}*beta{j};
    end
    e = e.*(rootprec*ones(1, size(e,2)));
    Sigma(1:2,1:2) = iwishrnd(e(:,1:2)'*e(:,1:2), N);

    % C. sampleGamma & sampleRho
    lmu = 3:(2+T); lop = 1:2;
    y_g = reshape(e(:,lmu), T*N, 1);
    x_g = repmat(e(:,lop), T, 1);
    ixx = (x_g'*x_g)\eye(size(x_g,2));
    [rchol, ~] = chol_c(ixx*gamma(3)*gamma(3)/(1-rho^2));
    gamma(1:2) = (randn(1,2)*rchol + (ixx*(x_g'*y_g))')';
    gamma(3) = 1/sqrt(gamrnd(N*T/2, 2/((1-rho^2)*sum((y_g - x_g*gamma(1:2)(:)).^2))));
    Sigma = combineSigOPbeta(Sigma(1:2,1:2), gamma, rho, T);

    % D. sampleSigL
    muL = u(:,3:end)';
    inc = data.choice > 0;
    ki = shape + sum(inc, 1)/2;
    thetai = zeros(size(ki));
    for i_obs = 1:N
      thetai(i_obs) = 2*theta./(2 + theta*sum(((loglambda(inc(:,i_obs),i_obs) - muL(inc(:,i_obs),i_obs)).*rootprec(i_obs)).^2));
    end
    sigL = sampleSigL_c(u(:,1), loglambda, ki, thetai, muL, sigL, exp(u(:,2)), ...
                        xint, wint, data.deduct, data.maxoop, data.prem, data.avail, ...
                        data.choice, lamlo, data.maxIter, data.threads, data.varLo, data.varHi, false);
    sigL = sigL(:);

    % E. sampleLamlo
    lam = exp(loglambda) - (lamlo*ones(1, T))';
    mull = data.xll*bll;
    [lamlo, loglambda] = sampleLamlo_c(lam, muL, sigL, mull, sigll, lamlo, data.choice, ...
                                       exp(u(:,1)), exp(u(:,2)), xint, wint, data.deduct, ...
                                       data.maxoop, data.prem, data.lamloHi, rootprec, data.threads, false);

    % Assert finite, valid state after each iteration
    assert(all(isfinite(B)), sprintf('Iter %d: Non-finite beta drawn', iter));
    assert(all(diag(Sigma) > 0), sprintf('Iter %d: Non-positive variance in Sigma', iter));
    assert(all(isfinite(sigL)), sprintf('Iter %d: Non-finite sigL drawn', iter));
    assert(all(isfinite(lamlo)), sprintf('Iter %d: Non-finite lamlo drawn', iter));
  end

  results.passed = results.passed + 1;
  fprintf('    [PASS] 3-iteration synthetic Gibbs cycle\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Synthetic Gibbs integration: %s', err.message);
  fprintf('    [FAIL] Synthetic Gibbs integration: %s\n', err.message);
end

end
