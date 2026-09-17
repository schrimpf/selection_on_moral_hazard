function results = test_sample_mu()
%% test_sample_mu.m
% Unit test for sampleMu_c conditional distribution calculations in panel data (T > 1).
% This test specifically checks the conditional Gaussian moments for period t > 1.
% It is expected to FAIL due to the indexing bug in sampleMu.c (lines 109, 113, 170)
% where period 1 indices (Sig[2*Tp2+2] and Tp2-1) are hardcoded instead of indexing period t.

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_sample_mu (checking conditional distribution for T=2)...\n');

[~, xint, wint] = gqzero(12);
xint = xint'*sqrt(2);
wint = wint/sqrt(pi);

%% 1. Test Period 1 (t=1): Should PASS because t=1 matches hardcoded indices
try
  T = 2;
  N = 10000;
  sigOP = [0.8 0.2; 0.2 0.5];
  gam = [0.3; 0.1; 0.8];
  rho = 0.7;
  Sig = combineSigOPbeta(sigOP, gam, rho, T);

  logomega = zeros(N, 1);
  logpsi   = zeros(N, 1);
  muL_init = [1.5 * ones(1, N); 0.5 * ones(1, N)];
  sigL     = 0.5 * ones(N, 1);
  lamlo    = zeros(N, 1);
  stdevi   = ones(N, 1);
  xbO      = zeros(N, 1);
  xbP      = zeros(N, 1);
  xbL      = zeros(T, N);
  loglambda = [1.0 * ones(1, N); 2.0 * ones(1, N)];

  np = 3;
  deduct = -ones(np, T, N);
  maxoop = -ones(np, T, N);
  prem   = zeros(np, T, N);
  avail  = zeros(np, T, N);
  choice = -ones(T, N); % bypass choice rejection to test raw Gaussian draws

  maxIter = 10;
  nthread = 4;
  multMH  = 0;

  muL_draws = sampleMu_c(logomega, loglambda, logpsi, muL_init, ...
                         xbO, xbP, xbL, sigL, Sig, ...
                         xint, wint, deduct, maxoop, prem, avail, choice, ...
                         lamlo, stdevi, maxIter, nthread, multMH);

  % Theoretical conditional mean for Period 1 (t=1):
  t = 1;
  idx_target = 2 + t;
  idx_u_other = [1, 2, setdiff(3:(2+T), idx_target)];
  V_theory = zeros(T+2, T+2);
  V_theory(1:(T+1), 1:(T+1)) = Sig(idx_u_other, idx_u_other);
  V_theory(1:(T+1), T+2)     = Sig(idx_u_other, idx_target);
  V_theory(T+2, 1:(T+1))     = Sig(idx_target, idx_u_other);
  V_theory(T+2, T+2)         = Sig(idx_target, idx_target) + sigL(1)^2;
  C_theory = [Sig(idx_target, idx_u_other)'; Sig(idx_target, idx_target)];
  var_theory_1 = Sig(idx_target, idx_target) - C_theory' * (V_theory \ C_theory);
  z_val_1 = [logomega(1); logpsi(1); muL_init(2,1); loglambda(1,1)];
  mean_theory_1 = xbL(1,1) + C_theory' * (V_theory \ z_val_1);

  emp_mean_1 = mean(muL_draws(1, :));
  emp_var_1  = var(muL_draws(1, :));
  se_mean_1  = sqrt(emp_var_1 / N);

  assert(~any(isnan(muL_draws(1, :))), 'Period 1 draws contain NaNs');
  assert(abs(emp_mean_1 - mean_theory_1) < 4 * se_mean_1, ...
         sprintf('Period 1 mean (%.4f) deviates from theory (%.4f)', emp_mean_1, mean_theory_1));

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleMu_c Period 1 (t=1) matches theoretical conditional Gaussian\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleMu_c Period 1: %s', err.message);
  fprintf('    [FAIL] sampleMu_c Period 1: %s\n', err.message);
end

%% 2. Test Period 2 (t=2) under AR(1) panel: verifies corrected conditional moments
try
  % Theoretical conditional moments for Period 2 (t=2):
  t = 2;
  idx_target = 2 + t;
  idx_u_other = [1, 2, setdiff(3:(2+T), idx_target)];
  V_theory = zeros(T+2, T+2);
  V_theory(1:(T+1), 1:(T+1)) = Sig(idx_u_other, idx_u_other);
  V_theory(1:(T+1), T+2)     = Sig(idx_u_other, idx_target);
  V_theory(T+2, 1:(T+1))     = Sig(idx_target, idx_u_other);
  V_theory(T+2, T+2)         = Sig(idx_target, idx_target) + sigL(1)^2;
  C_theory = [Sig(idx_target, idx_u_other)'; Sig(idx_target, idx_target)];
  var_theory_2 = Sig(idx_target, idx_target) - C_theory' * (V_theory \ C_theory);
  weights = C_theory' / V_theory;

  % In Gibbs sampling, Period 2 conditions on the newly updated muL_draws(1, :)
  mean_theory_2 = weights * [0; 0; mean(muL_draws(1, :)); 2.0];
  expected_var_2 = var_theory_2 + weights(3)^2 * var(muL_draws(1, :));

  emp_mean_2 = mean(muL_draws(2, :));
  emp_var_2  = var(muL_draws(2, :));
  se_mean_2  = sqrt(emp_var_2 / N);

  assert(~any(isnan(muL_draws(2, :))), 'Period 2 draws contain NaNs');
  assert(abs(emp_mean_2 - mean_theory_2) < 4 * se_mean_2, ...
         sprintf('Period 2 mean (%.4f) deviates from theory (%.4f)', emp_mean_2, mean_theory_2));
  assert(abs(emp_var_2 - expected_var_2) < 0.05, ...
         sprintf('Period 2 variance (%.4f) deviates from theory (%.4f)', emp_var_2, expected_var_2));

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleMu_c Period 2 (t=2) matches theoretical conditional Gaussian\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleMu_c Period 2 AR(1): %s', err.message);
  fprintf('    [FAIL] sampleMu_c Period 2 AR(1): %s\n', err.message);
end

%% 3. Test Period 2 (t=2) with unequal period variances: verifies stability and correct moments
try
  % Covariance matrix with Var(muL_1) = 1.0 and Var(muL_2) = 4.0
  Sig_unequal = [ 1.0   0.2   0.3   0.4;
                  0.2   1.0   0.1   0.2;
                  0.3   0.1   1.0   1.0;
                  0.4   0.2   1.0   4.0 ];

  muL_draws_unequal = sampleMu_c(logomega, loglambda, logpsi, muL_init, ...
                                 xbO, xbP, xbL, sigL, Sig_unequal, ...
                                 xint, wint, deduct, maxoop, prem, avail, choice, ...
                                 lamlo, stdevi, maxIter, nthread, multMH);

  assert(~any(isnan(muL_draws_unequal(1, :))), 'Period 1 draws contain NaNs');
  assert(~any(isnan(muL_draws_unequal(2, :))), 'Period 2 draws contain NaNs');

  t = 2;
  idx_target = 2 + t;
  idx_u_other = [1, 2, setdiff(3:(2+T), idx_target)];
  V_u = zeros(T+2, T+2);
  V_u(1:(T+1), 1:(T+1)) = Sig_unequal(idx_u_other, idx_u_other);
  V_u(1:(T+1), T+2)     = Sig_unequal(idx_u_other, idx_target);
  V_u(T+2, 1:(T+1))     = Sig_unequal(idx_target, idx_u_other);
  V_u(T+2, T+2)         = Sig_unequal(idx_target, idx_target) + sigL(1)^2;
  C_u = [Sig_unequal(idx_target, idx_u_other)'; Sig_unequal(idx_target, idx_target)];
  var_theory_u2 = Sig_unequal(idx_target, idx_target) - C_u' * (V_u \ C_u);
  weights_u = C_u' / V_u;

  mean_theory_u2 = weights_u * [0; 0; mean(muL_draws_unequal(1, :)); 2.0];
  expected_var_u2 = var_theory_u2 + weights_u(3)^2 * var(muL_draws_unequal(1, :));

  emp_mean_u2 = mean(muL_draws_unequal(2, :));
  emp_var_u2  = var(muL_draws_unequal(2, :));
  se_mean_u2  = sqrt(emp_var_u2 / N);

  assert(abs(emp_mean_u2 - mean_theory_u2) < 4 * se_mean_u2, ...
         sprintf('Unequal Var Period 2 mean (%.4f) deviates from theory (%.4f)', emp_mean_u2, mean_theory_u2));
  assert(abs(emp_var_u2 - expected_var_u2) < 0.05, ...
         sprintf('Unequal Var Period 2 variance (%.4f) deviates from theory (%.4f)', emp_var_u2, expected_var_u2));

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleMu_c Period 2 (t=2) unequal variances matches theory\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleMu_c Period 2 unequal variances: %s', err.message);
  fprintf('    [FAIL] sampleMu_c Period 2 unequal variances: %s\n', err.message);
end

end
