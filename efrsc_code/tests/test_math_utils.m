function results = test_math_utils()
%% test_math_utils.m
% Unit tests for mathematical utilities and distribution samplers.

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_math_utils...\n');

%% 1. Test chol_c (Cholesky decomposition MEX)
try
  % Positive definite test
  A = [4 2 1; 2 5 3; 1 3 6];
  [R, fail] = chol_c(A);
  assert(fail == 0, 'chol_c reported failure on positive-definite matrix');
  assert(max(abs(R'*R - A)(:)) < 1e-12, 'chol_c R''*R does not reconstruct A');
  assert(all(diag(R) > 0), 'chol_c diagonal elements must be positive');

  % Identity matrix test
  I5 = eye(5);
  [R_eye, fail_eye] = chol_c(I5);
  assert(fail_eye == 0, 'chol_c failed on identity');
  assert(max(abs(R_eye - I5)(:)) < 1e-14, 'chol_c did not return identity for I');

  % Non-positive definite matrix test
  A_bad = [1 2; 2 1]; % eigenvalues: 3, -1
  [~, fail_bad] = chol_c(A_bad);
  assert(fail_bad ~= 0, 'chol_c did not detect indefinite matrix');

  results.passed = results.passed + 1;
  fprintf('    [PASS] chol_c\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('chol_c: %s', err.message);
  fprintf('    [FAIL] chol_c: %s\n', err.message);
end

%% 2. Test randdtn_c (Truncated Normal Sampler)
try
  n_draws = 10000;
  mu = 1.0;
  sigma = 0.5; % standard deviation (randdtn_c squares it internally)
  a = 0.5;
  b = 2.0;

  % Vectorized draw test
  draws = randdtn_c(mu * ones(1, n_draws), sigma * ones(1, n_draws), ...
                    a * ones(1, n_draws), b * ones(1, n_draws));

  % Verify all draws strictly within bounds [a, b]
  assert(all(draws >= a), 'randdtn_c draw below lower bound');
  assert(all(draws <= b), 'randdtn_c draw above upper bound');

  % Analytical mean of truncated normal
  alpha = (a - mu) / sigma;
  beta = (b - mu) / sigma;
  Z = normcdf(beta) - normcdf(alpha);
  expected_mean = mu + (normpdf(alpha) - normpdf(beta)) / Z * sigma;
  sample_mean = mean(draws);

  % Check within 3.5 standard errors
  se = std(draws) / sqrt(n_draws);
  assert(abs(sample_mean - expected_mean) < 3.5 * se, ...
         sprintf('randdtn_c mean (%.4f) deviates from analytical mean (%.4f)', ...
                 sample_mean, expected_mean));

  results.passed = results.passed + 1;
  fprintf('    [PASS] randdtn_c\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('randdtn_c: %s', err.message);
  fprintf('    [FAIL] randdtn_c: %s\n', err.message);
end

%% 3. Test gqzero (Gauss-Hermite Quadrature)
try
  order = 24;
  [ord, xint, wint] = gqzero(order);
  assert(ord == order, 'gqzero returned incorrect order');
  assert(numel(xint) == numel(wint), 'gqzero nodes and weights size mismatch');
  assert(all(wint > 0), 'gqzero weights must be positive');

  % Check symmetry of nodes and weights
  assert(max(abs(xint + flip(xint))) < 1e-12, 'gqzero nodes must be symmetric around 0');
  assert(max(abs(wint - flip(wint))) < 1e-12, 'gqzero weights must be symmetric');

  % Integral of exp(-x^2) = sqrt(pi) -> sum(wint) = sqrt(pi)
  assert(abs(sum(wint) - sqrt(pi)) < 1e-12, 'Quadrature 0-th moment must equal sqrt(pi)');

  % Integral of x^2 * exp(-x^2) = sqrt(pi) / 2
  assert(abs(sum(wint .* (xint.^2)) - sqrt(pi)/2) < 1e-11, ...
         'Quadrature 2-nd moment must equal sqrt(pi)/2');

  results.passed = results.passed + 1;
  fprintf('    [PASS] gqzero\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('gqzero: %s', err.message);
  fprintf('    [FAIL] gqzero: %s\n', err.message);
end

%% 4. Test iwishrnd (Inverse Wishart RNG)
try
  p = 3;
  df = 20;
  Psi = [2 0.5 0.2; 0.5 1.5 0.3; 0.2 0.3 1.0];

  % Single draw test
  W = iwishrnd(Psi, df);
  assert(all(size(W) == [p p]), 'iwishrnd output dimension mismatch');
  assert(max(abs(W - W')) < 1e-12, 'iwishrnd output must be symmetric');
  [~, p_fail] = chol_c(W);
  assert(p_fail == 0, 'iwishrnd output must be positive definite');

  % Mean of Inverse Wishart is Psi / (df - p - 1)
  n_mc = 2000;
  sum_W = zeros(p, p);
  for k = 1:n_mc
    sum_W = sum_W + iwishrnd(Psi, df);
  end
  mean_W = sum_W / n_mc;
  expected_W = Psi / (df - p - 1);
  diff_norm = norm(mean_W - expected_W, 'fro') / norm(expected_W, 'fro');
  assert(diff_norm < 0.10, sprintf('iwishrnd mean relative error (%.3f) too high', diff_norm));

  results.passed = results.passed + 1;
  fprintf('    [PASS] iwishrnd\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('iwishrnd: %s', err.message);
  fprintf('    [FAIL] iwishrnd: %s\n', err.message);
end

%% 5. Test combineSigOPbeta
try
  SigmaOP = [0.8 0.2; 0.2 0.6];
  gamma = [0.3; 0.1; 0.5]; % gamma(1:2) regression, gamma(3) std dev
  rho = 0.4;
  T = 3;

  Sig = combineSigOPbeta(SigmaOP, gamma, rho, T);
  expected_dim = 2 + T;
  assert(all(size(Sig) == [expected_dim expected_dim]), 'combineSigOPbeta size mismatch');
  assert(max(abs(Sig - Sig')) < 1e-12, 'combineSigOPbeta must produce symmetric matrix');
  [~, p_fail] = chol_c(Sig);
  assert(p_fail == 0, 'combineSigOPbeta must produce positive definite matrix');
  assert(max(abs(Sig(1:2, 1:2) - SigmaOP)(:)) < 1e-14, 'OP block must match SigmaOP');

  results.passed = results.passed + 1;
  fprintf('    [PASS] combineSigOPbeta\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('combineSigOPbeta: %s', err.message);
  fprintf('    [FAIL] combineSigOPbeta: %s\n', err.message);
end

%% 6. Test sampleShape_c
try
  shape0 = 2.5;
  theta0 = 1.2;
  n_obs = 100;
  x = gamrnd(shape0, theta0, n_obs, 1);
  varHi = 100.0;

  shape_new = sampleShape_c(shape0, theta0, x, 1/varHi);
  assert(isfinite(shape_new), 'sampleShape_c returned non-finite value');
  assert(shape_new > 0, 'sampleShape_c returned non-positive value');

  results.passed = results.passed + 1;
  fprintf('    [PASS] sampleShape_c\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('sampleShape_c: %s', err.message);
  fprintf('    [FAIL] sampleShape_c: %s\n', err.message);
end

end
