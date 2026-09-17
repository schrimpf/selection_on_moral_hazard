function results = test_latent_samplers()
%% test_latent_samplers.m
% Unit tests for conditional latent variable samplers.

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_latent_samplers...\n');

[~, xint, wint] = gqzero(24);
xint = xint'*sqrt(2);
wint = wint/sqrt(pi);

%% 1. Test sampleSigL_c (Idiosyncratic Variance Sampler)
try
  N = 30;
  T = 2;
  np = 3;

  deduct = 500 * rand(np, T, N);
  maxoop = 1500 * rand(np, T, N) + 500;
  prem = 100 * rand(np, T, N);
  avail = ones(np, T, N);

  logomega = zeros(N, 1);
  psi = 0.001 * ones(N, 1);
  muL = 6.0 * ones(T, N);
  sigL = 0.5 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = muL + 0.1 * randn(T, N);
  totalSpend = zeros(T, N);

  % Generate valid choices for baseline
  [choice, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                              sigL, psi, xint, wint, deduct, maxoop, ...
                              prem, avail, ones(T, N), totalSpend, 100);

  shape = 2.0;
  theta = 1.0;
  inc = choice > 0;
  ki = shape + sum(inc, 1) / 2;
  thetai = zeros(size(ki));
  for i = 1:N
    thetai(i) = 2*theta / (2 + theta*sum((loglambda(inc(:,i), i) - muL(inc(:,i), i)).^2));
  end

  varLo = 0.01;
  varHi = 10.0;
  maxIter = 200;
  threads = 4;

  sigL_new = sampleSigL_c(logomega, loglambda, ki, thetai, muL, sigL, psi, ...
                          xint, wint, deduct, maxoop, prem, avail, ...
                          choice, lamlo, maxIter, threads, varLo, varHi, false);

  assert(numel(sigL_new) == N, 'sampleSigL_c returned incorrect number of elements');
  assert(all(isfinite(sigL_new)), 'sampleSigL_c returned non-finite values');
  assert(all(sigL_new(:) >= sqrt(varLo) & sigL_new(:) <= sqrt(varHi)), ...
         'sampleSigL_c returned draws outside variance bounds [varLo, varHi]');

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleSigL_c (support bounds & finite values)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleSigL_c: %s', err.message);
  fprintf('    [FAIL] sampleSigL_c: %s\n', err.message);
end

%% 2. Test samplePsi_c (Risk Aversion Sampler)
try
  N = 30;
  T = 2;
  np = 3;

  deduct = 500 * rand(np, T, N);
  maxoop = 1500 * rand(np, T, N) + 500;
  prem = 100 * rand(np, T, N);
  avail = ones(np, T, N);

  logomega = zeros(N, 1);
  psi = 0.001 * ones(N, 1);
  logpsi = log(psi);
  muL = 6.0 * ones(T, N);
  sigL = 0.5 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = muL;
  totalSpend = zeros(T, N);

  [choice, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                              sigL, psi, xint, wint, deduct, maxoop, ...
                              prem, avail, ones(T, N), totalSpend, 100);

  muP = logpsi;
  sigP = 0.3 * ones(N, 1);
  stdevi = ones(N, 1);
  maxIter = 200;
  threads = 4;

  logpsi_new = samplePsi_c(logomega, loglambda, muP, sigP, muL, ...
                           sigL, logpsi, xint, wint, deduct, maxoop, prem, ...
                           avail, choice, lamlo, stdevi, maxIter, threads, false);

  assert(numel(logpsi_new) == N, 'samplePsi_c dimension mismatch');
  assert(all(isfinite(logpsi_new)), 'samplePsi_c returned non-finite values');
  psi_new = exp(logpsi_new);
  assert(all(psi_new > 0), 'samplePsi_c risk aversion must be strictly positive');

  results.passed = results.passed + 1;
  fprintf('    [PASS] samplePsi_c (positivity & finite draws)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('samplePsi_c: %s', err.message);
  fprintf('    [FAIL] samplePsi_c: %s\n', err.message);
end

%% 3. Test sampleLamlo_c (Minimum Spending Bound Sampler)
try
  N = 30;
  T = 2;
  np = 3;

  deduct = 500 * rand(np, T, N);
  maxoop = 1500 * rand(np, T, N) + 500;
  prem = 100 * rand(np, T, N);
  avail = ones(np, T, N);

  omega = ones(N, 1);
  psi = 0.001 * ones(N, 1);
  muL = 6.0 * ones(T, N);
  sigL = 0.5 * ones(N, 1);
  lamlo0 = 50 * ones(N, 1);
  lam = exp(muL);
  loglambda = log(lam + (lamlo0 * ones(1, T))');
  totalSpend = zeros(T, N);

  [choice, ~] = findChoices_c(log(omega), loglambda, lamlo0, 1, muL, ...
                              sigL, psi, xint, wint, deduct, maxoop, ...
                              prem, avail, ones(T, N), totalSpend, 100);

  mull = 50 * ones(N, 1);
  sigll = 10.0;
  lamloHi = 2000.0;
  rootprec = ones(N, 1);
  threads = 4;

  [lamlo_new, loglambda_new] = sampleLamlo_c(lam, muL, sigL, mull, sigll, lamlo0, choice, ...
                                             omega, psi, xint, wint, deduct, ...
                                             maxoop, prem, lamloHi, rootprec, threads, false);

  assert(numel(lamlo_new) == N, 'sampleLamlo_c dimension mismatch');
  assert(all(isfinite(lamlo_new)), 'sampleLamlo_c returned non-finite values');
  assert(all(lamlo_new <= lamloHi), 'sampleLamlo_c exceeded lamloHi');
  assert(all(isfinite(loglambda_new(:))), 'loglambda_new has non-finite values');

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleLamlo_c (upper bound & finite draws)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleLamlo_c: %s', err.message);
  fprintf('    [FAIL] sampleLamlo_c: %s\n', err.message);
end

%% 4. Test findValidLatent3_c with badindex targeting
try
  N = 10;
  T = 2;
  np = 3;

  deduct = 500 * rand(np, T, N);
  maxoop = 1500 * rand(np, T, N) + 500;
  prem = 100 * rand(np, T, N);
  avail = ones(np, T, N);

  logomega = zeros(N, 1);
  psi = 0.001 * ones(N, 1);
  muL = 6.0 * ones(T, N);
  sigL = 0.5 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = muL;
  totalSpend = zeros(T, N);

  % Baseline choices
  [choice, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                              sigL, psi, xint, wint, deduct, maxoop, ...
                              prem, avail, ones(T, N), totalSpend, 100);

  % Bounds setup
  NU = 4 + T;
  lb = -10 * ones(NU, N);
  ub =  10 * ones(NU, N);
  lb(1, :) = log(100); ub(1, :) = log(50000);
  lb(2, :) = log(1e-6); ub(2, :) = log(0.05);
  lb(3:2+T, :) = 2; ub(3:2+T, :) = 12;
  lb(3+T, :) = log(0.01); ub(3+T, :) = log(10);
  lb(4+T, :) = 0; ub(4+T, :) = 2000;

  % Intentionally corrupt 2 observations (indices 3 and 7)
  bad_target = [3, 7];
  muL_test = muL;
  muL_test(:, bad_target) = -50.0; % extreme value violating choices

  % Call findValidLatent3_c with badindex
  threads = 4;
  [u_out, sigL_out, lamlo_out, fail] = findValidLatent3_c( ...
      muL_test, log(sigL), logomega, log(psi), ...
      deduct, maxoop, prem, xint, wint, choice, threads, ...
      lb, ub, lamlo, totalSpend, bad_target);

  % Check output shapes
  assert(all(size(u_out) == [2+T, N]), 'u_out shape mismatch');
  assert(numel(sigL_out) == N, 'sigL_out length mismatch');
  assert(numel(lamlo_out) == N, 'lamlo_out length mismatch');

  % Non-targeted observations (e.g. index 1) must be preserved exactly
  assert(abs(u_out(1, 1) - logomega(1)) < 1e-12, 'Non-bad index logomega modified');
  assert(abs(sigL_out(1) - sigL(1)) < 1e-12, 'Non-bad index sigL modified');

  results.passed = results.passed + 1;
  fprintf('    [PASS] findValidLatent3_c (selective badindex targeting)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('findValidLatent3_c: %s', err.message);
  fprintf('    [FAIL] findValidLatent3_c: %s\n', err.message);
end

end
