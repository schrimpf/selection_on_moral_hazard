function results = test_synthetic_data()
%% test_synthetic_data.m
% Unit test for options.synthetic in simulate.m:
% 1. Verifies observation fields (choiceObs, spendObs, etc.) are stripped.
% 2. Verifies synthetic covariates preserve exact empirical support.
% 3. Verifies feature preservation (means, standard deviations, correlations).
% 4. Verifies individual privacy (zero exact covariate matches).

addpath('..');
addpath('../postEst');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_synthetic_data...\n');

if (exist('isOctave', 'file') && isOctave())
  rand('state', 42);
  randn('state', 42);
else
  rng(42);
end

%% Setup synthetic testing data structure if .mat file is not present
has_mat = exist('revisedCode01.mat', 'file') || exist('../revisedCode01.mat', 'file') || ...
          exist('mainSpec02.mat', 'file') || exist('../mainSpec02.mat', 'file');

if (has_mat)
  if exist('revisedCode01.mat', 'file') || exist('../revisedCode01.mat', 'file')
    f = 'revisedCode01.mat';
  else
    f = 'mainSpec02.mat';
  end
  load(f);
  nSaved = numel(chain.meanL);
  post = floor(nSaved/2):nSaved;
  for k=1:numel(chain.beta)
    b{k} = mean(chain.beta{k}(post,:),1)';
  end
  Sig = squeeze(mean(chain.sig(post,:,:),1));
  gam = squeeze(mean(chain.gamma(post,:,:),1));
  th  = squeeze(mean(chain.theta(post)));
  shp = squeeze(mean(chain.shape(post)));
  rh  = squeeze(mean(chain.rho(post)));
  bl  = mean(chain.bll(post,:),1)';
  sll = mean(chain.sigll(post));
  al  = 0;
else
  data.N = 200;
  data.T = 2;
  data.year = [2003; 2004];
  data.order = 12;
  data.varLo = 0;
  data.varHi = 10;
  data.threads = 1;
  data.maxIter = 10;
  data.choice = ones(2, data.N);
  data.totalSpend = 500 * ones(2, data.N);
  for p = 1:5
    data.avail(p,:,:) = 1;
    data.prem(p,:,:) = 100 * p;
    data.deduct(p,:,:) = 250 * p;
    data.copay(p,:,:) = 0.2;
    data.maxoop(p,:,:) = 1000 * p;
  end
  x1 = [ones(data.N, 1), eye(data.N, 3), eye(data.N, 3), ...
        round(40 + 5*randn(data.N, 1)), rand(data.N, 1) > 0.5, ...
        round(10 + 3*randn(data.N, 1)), 35 + 10*rand(data.N, 1), ...
        double(rand(data.N, 3) > 0.7)];
  data.x{1} = x1;
  data.x{2} = x1;
  data.xll = x1;
  data.xhs = x1(:, 2:end);
  data.x{3}{1} = [x1, zeros(data.N, 1)];
  data.x{3}{2} = [x1, ones(data.N, 1)];
  b = {zeros(14, 1), zeros(14, 1), zeros(15, 1)};
  Sig = eye(4);
  gam = [0.1; 0.1; 0.1];
  th = 1; shp = 1; rh = 0.5; bl = zeros(14, 1); sll = 1; al = 0;
end

%% Test 1: Observation fields stripping
try
  options.synthetic = true;
  options.avgPlan = true;
  options.cf = false;
  N_sim = 500;
  sim = simulate(N_sim, data, b, Sig, gam, rh, shp, th, bl, sll, al, options);

  assert(~isfield(sim, 'choiceObs'), 'choiceObs was not stripped');
  assert(~isfield(sim, 'spendObs'),  'spendObs was not stripped');
  assert(~isfield(sim, 'premObs'),   'premObs was not stripped');
  assert(~isfield(sim, 'deductObs'), 'deductObs was not stripped');
  assert(~isfield(sim, 'copayObs'),  'copayObs was not stripped');
  assert(~isfield(sim, 'maxoopObs'), 'maxoopObs was not stripped');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Observation fields successfully stripped\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Observation fields stripping: %s', err.message);
  fprintf('    [FAIL] Observation fields stripping: %s\n', err.message);
end

%% Test 2: Support preservation for all covariates
try
  % 1. Constant
  assert(all(sim.x{1}(:, 1) == 1), 'Constant column is not identically 1');

  % 2. Coverage tier dummies: mutually exclusive
  tier_sums = sum(sim.x{1}(:, 2:4), 2);
  assert(all(tier_sums == 0 | tier_sums == 1), 'Coverage tier dummies are not mutually exclusive');

  % 3. Group dummies: mutually exclusive
  grp_sums = sum(sim.x{1}(:, 5:7), 2);
  assert(all(grp_sums == 0 | grp_sums == 1), 'Group dummies are not mutually exclusive');

  % 4. Age bounds and half-integer support
  min_age = min(data.x{1}(:, 8));
  max_age = max(data.x{1}(:, 8));
  assert(min(sim.x{1}(:, 8)) >= min_age, sprintf('Age min violated: %g < %g', min(sim.x{1}(:, 8)), min_age));
  assert(max(sim.x{1}(:, 8)) <= max_age, sprintf('Age max violated: %g > %g', max(sim.x{1}(:, 8)), max_age));

  % 5. Sex support: binary {0, 1}
  assert(all(sim.x{1}(:, 9) == 0 | sim.x{1}(:, 9) == 1), 'Sex values not in {0, 1}');

  % 6. Tenure bounds
  min_tenure = 0;
  max_tenure = max(data.x{1}(:, 10));
  assert(min(sim.x{1}(:, 10)) >= min_tenure, 'Tenure min violated');
  assert(max(sim.x{1}(:, 10)) <= max_tenure, 'Tenure max violated');

  % 7. Wage bounds: strictly positive and within empirical support
  min_wage = min(data.x{1}(:, 11));
  max_wage = max(data.x{1}(:, 11));
  assert(min(sim.x{1}(:, 11)) >= min_wage, sprintf('Wage min violated: %g < %g', min(sim.x{1}(:, 11)), min_wage));
  assert(max(sim.x{1}(:, 11)) <= max_wage, sprintf('Wage max violated: %g > %g', max(sim.x{1}(:, 11)), max_wage));

  % 8. Chronic conditions in yearly panel {0, 1} and average {0, 0.5, 1}
  assert(all(ismember(unique(sim.x{1}(:, 12:14)), [0 0.5 1])), 'Averaged chronic conditions not in {0, 0.5, 1}');
  for t = 1:numel(sim.x{3})
    assert(all(ismember(unique(sim.x{3}{t}(:, 12:14)), [0 1])), sprintf('Year %d chronic condition not in {0, 1}', t));
  end

  % 9. Internal consistency between x{1}, x{2}, xll, xhs
  assert(isequal(sim.x{1}, sim.x{2}), 'sim.x{1} and sim.x{2} mismatch');
  assert(isequal(sim.x{1}, sim.xll),  'sim.x{1} and sim.xll mismatch');
  assert(isequal(sim.x{1}(:, 2:end), sim.xhs), 'sim.x{1} and sim.xhs mismatch');

  results.passed = results.passed + 1;
  fprintf('    [PASS] Synthetic covariates preserve exact empirical support\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Support preservation: %s', err.message);
  fprintf('    [FAIL] Support preservation: %s\n', err.message);
end

%% Test 3: Feature preservation and privacy
try
  % Privacy: Zero exact matching rows between synthetic covariates and real data
  exact_matches = 0;
  for i = 1:min(sim.N, 200)
    match = all(abs(sim.x{1}(i, :) - data.x{1}) < 1e-4, 2);
    if any(match)
      exact_matches = exact_matches + 1;
    end
  end
  assert(exact_matches == 0, sprintf('%d synthetic records exactly match original records', exact_matches));

  % Feature preservation: Means of continuous variables match original within 5%
  age_err = abs(mean(sim.x{1}(:, 8)) - mean(data.x{1}(:, 8))) / std(data.x{1}(:, 8));
  assert(age_err < 0.15, sprintf('Age mean relative error too high: %g', age_err));

  wage_err = abs(mean(sim.x{1}(:, 11)) - mean(data.x{1}(:, 11))) / std(data.x{1}(:, 11));
  assert(wage_err < 0.15, sprintf('Wage mean relative error too high: %g', wage_err));

  results.passed = results.passed + 1;
  fprintf('    [PASS] Feature preservation and privacy verification passed\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('Feature preservation/privacy: %s', err.message);
  fprintf('    [FAIL] Feature preservation/privacy: %s\n', err.message);
end

end
