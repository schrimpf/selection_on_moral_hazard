function results = test_choice_engine()
%% test_choice_engine.m
% Unit tests for choice evaluation engine (findChoices_c, choiceFnEx).

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_choice_engine...\n');

%% Setup quadrature and common parameters
[~, xint, wint] = gqzero(24);
xint = xint'*sqrt(2);
wint = wint/sqrt(pi);

%% 1. Monotonicity test: Strictly better plan (lower deductible & max OOP) must be preferred
try
  N = 10;
  T = 1;
  np = 3;

  % 3 plans:
  % Plan 1: High deductible ($1000), High max OOP ($3000), Prem = 0
  % Plan 2: Medium deductible ($500), Medium max OOP ($1500), Prem = 0
  % Plan 3: Low deductible ($100), Low max OOP ($500), Prem = 0
  % With zero premiums across all plans, Plan 3 strictly dominates Plan 2, which dominates Plan 1!
  deduct = zeros(np, T, N);
  maxoop = zeros(np, T, N);
  prem = zeros(np, T, N);
  avail = ones(np, T, N);

  for i = 1:N
    deduct(:, 1, i) = [1000; 500; 100];
    maxoop(:, 1, i) = [3000; 1500; 500];
  end

  logomega = zeros(N, 1);       % omega = 1
  psi = 0.001 * ones(N, 1);      % risk aversion
  logpsi = log(psi);
  muL = 6.0 * ones(T, N);        % E[lambda] around 400
  sigL = 0.5 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = (randn(N, T).*(sigL*ones(1, T)) + muL')';
  totalSpend = zeros(T, N);
  choice_init = ones(T, N);

  [choices, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                               sigL, psi, xint, wint, deduct, maxoop, ...
                               prem, avail, choice_init, totalSpend, 100);

  % All individuals must strictly choose Plan 3 (1-indexed choice == 3)
  assert(all(choices == 3), 'All individuals must choose strictly dominant plan (Plan 3)');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Dominant plan selection (monotonicity)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Monotonicity test: %s', err.message);
  fprintf('    [FAIL] Dominant plan selection: %s\n', err.message);
end

%% 2. Premium elasticity test: High premium forces individual into lower-cost plan
try
  N = 10;
  T = 1;
  np = 2;

  % Plan 1: Cheap plan (Deduct = 1000, MaxOOP = 2000, Prem = 0)
  % Plan 2: Expensive generous plan (Deduct = 100, MaxOOP = 500, Prem = 50,000)
  % With an exorbitant premium on Plan 2, all individuals must choose Plan 1!
  deduct = zeros(np, T, N);
  maxoop = zeros(np, T, N);
  prem = zeros(np, T, N);
  avail = ones(np, T, N);

  for i = 1:N
    deduct(:, 1, i) = [1000; 100];
    maxoop(:, 1, i) = [2000; 500];
    prem(:, 1, i)   = [0; 50000];
  end

  logomega = zeros(N, 1);
  psi = 0.001 * ones(N, 1);
  muL = 5.0 * ones(T, N);
  sigL = 0.5 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = muL;
  totalSpend = zeros(T, N);
  choice_init = ones(T, N);

  [choices, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                               sigL, psi, xint, wint, deduct, maxoop, ...
                               prem, avail, choice_init, totalSpend, 100);

  assert(all(choices == 1), 'Exorbitant premium must cause individuals to reject Plan 2');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Premium elasticity\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Premium elasticity: %s', err.message);
  fprintf('    [FAIL] Premium elasticity: %s\n', err.message);
end

%% 3. Risk aversion sensitivity: Higher psi shifts preference toward more comprehensive insurance
try
  np = 2;
  T = 1;
  N = 2;

  % Plan 1: High deductible / low premium (Deduct = 2000, MaxOOP = 4000, Prem = 0)
  % Plan 2: Low deductible / modest premium (Deduct = 0, MaxOOP = 500, Prem = 1500)
  deduct = zeros(np, T, N);
  maxoop = zeros(np, T, N);
  prem = zeros(np, T, N);
  avail = ones(np, T, N);

  for i = 1:N
    deduct(:, 1, i) = [2000; 0];
    maxoop(:, 1, i) = [4000; 500];
    prem(:, 1, i)   = [0; 1500];
  end

  logomega = zeros(N, 1);
  % Individual 1: Risk neutral (psi = 1e-6) -> prefers cheaper high-deductible plan (Plan 1)
  % Individual 2: Highly risk averse (psi = 0.02) -> prefers comprehensive plan (Plan 2)
  psi = [1e-6; 0.02];
  muL = 7.0 * ones(T, N); % E[lambda] ~ 1100
  sigL = 0.8 * ones(N, 1); % high variance in spending risk
  lamlo = zeros(N, 1);
  loglambda = muL;
  totalSpend = zeros(T, N);
  choice_init = ones(T, N);

  [choices, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                               sigL, psi, xint, wint, deduct, maxoop, ...
                               prem, avail, choice_init, totalSpend, 100);

  assert(choices(1) == 1, 'Risk-neutral individual should prefer low-premium plan (Plan 1)');
  assert(choices(2) == 2, 'Highly risk-averse individual should prefer comprehensive coverage (Plan 2)');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Risk aversion sensitivity\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Risk aversion sensitivity: %s', err.message);
  fprintf('    [FAIL] Risk aversion sensitivity: %s\n', err.message);
end

%% 4. Multi-period panel consistency
try
  N = 20;
  T = 3;
  np = 4;

  deduct = 500 * rand(np, T, N);
  maxoop = 1500 * rand(np, T, N) + 500;
  prem = 200 * rand(np, T, N);
  avail = ones(np, T, N);

  logomega = 0.1 * randn(N, 1);
  psi = 0.001 * ones(N, 1);
  muL = 6.0 * ones(T, N);
  sigL = 0.4 * ones(N, 1);
  lamlo = zeros(N, 1);
  loglambda = muL;
  totalSpend = zeros(T, N);
  choice_init = ones(T, N);

  [choices, ~] = findChoices_c(logomega, loglambda, lamlo, 1, muL, ...
                               sigL, psi, xint, wint, deduct, maxoop, ...
                               prem, avail, choice_init, totalSpend, 100);

  assert(all(size(choices) == [T N]), 'Output choices dimensions mismatch [T, N]');
  assert(all(choices(:) >= 1 & choices(:) <= np), 'Choices must be valid plan indices 1..np');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Multi-period panel dimensions\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Multi-period test: %s', err.message);
  fprintf('    [FAIL] Multi-period test: %s\n', err.message);
end

%% 5. findChoicesNMH_c basic functionality and N=1 handling
try
  N = 1;
  T = 1;
  np = 3;
  deduct = [250; 500; 1000];
  maxoop = [1000; 2000; 3000];
  prem   = [1200; 800; 400];
  avail  = [1; 1; 1];
  choice = 0;
  spend  = 0;
  [c_nmh, v_nmh] = findChoicesNMH_c(0, 0, 0, 0.5, 0, 0.5, 1e-4, ...
                                    xint, wint, deduct, maxoop, prem, avail, choice, spend, 10);
  assert(c_nmh == 3, sprintf('NMH choice for N=1 expected 3, got %d', c_nmh));
  assert(size(v_nmh, 1) == np, 'NMH value dimension mismatch');

  results.passed = results.passed + 1;
  fprintf('    [PASS] findChoicesNMH_c N=1 singleton dimension support\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('findChoicesNMH N=1: %s', err.message);
  fprintf('    [FAIL] findChoicesNMH N=1: %s\n', err.message);
end

%% 6. findChoicesNMH_c MPFR fallback with excluded plan (deduct(1) = -1)
try
  N = 4;
  T = 1;
  np = 3;
  deduct = repmat([-1; 500; 1000], [1, T, N]);
  maxoop = repmat([1000; 2000; 3000], [1, T, N]);
  prem   = repmat([1200; 800; 400], [1, T, N]);
  avail  = repmat([0; 1; 1], [1, T, N]);
  choice = zeros(T, N);
  spend  = zeros(T, N);
  % High psi triggers overflow into choiceFnNMHuge
  psi = 1000 * ones(N, 1);
  [c_huge, v_huge] = findChoicesNMH_c(zeros(N,1), zeros(N,1), zeros(N,1), 0.5, ...
                                      zeros(T,N), 0.5*ones(N,1), psi, ...
                                      xint, wint, deduct, maxoop, prem, avail, choice, spend, ...
                                      10, 0, 1);
  assert(all(c_huge == 3), 'NMHuge failed to select valid plan when plan 1 is excluded');
  assert(all(v_huge(1,:,:) == -Inf), 'NMHuge expected -Inf for unavailable plan 1');

  results.passed = results.passed + 1;
  fprintf('    [PASS] findChoicesNMH_c MPFR with excluded plan & 18 args\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('findChoicesNMH MPFR: %s', err.message);
  fprintf('    [FAIL] findChoicesNMH MPFR: %s\n', err.message);
end

end
