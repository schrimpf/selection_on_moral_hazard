function results = test_matrix_kernels()
%% test_matrix_kernels.m
% Unit tests for matrix operations and MEX kernels (XtTimesKronSIdTimesY).

addpath('..');
results.passed = 0;
results.failed = 0;
results.errors = {};

fprintf('  Running test_matrix_kernels...\n');

%% 1. Standard dimension test
try
  N = 100;
  Ns = 4;
  Nx = 25;
  Ny = 25;

  x = randn(N*Ns, Nx);
  y = randn(N*Ns, Ny);
  S = randn(Ns, Ns);
  S = S + S'; % symmetric S

  % Reference via block linear combinations
  Z = zeros(N*Ns, Ny);
  for a = 1:Ns
    for b = 1:Ns
      Z((a-1)*N+1:a*N, :) = Z((a-1)*N+1:a*N, :) + S(a,b) * y((b-1)*N+1:b*N, :);
    end
  end
  prod_ref = x' * Z;

  % MEX kernel
  prod_mex = XtTimesKronSIdTimesY(x, S, N, y);

  max_diff = max(abs(prod_mex(:) - prod_ref(:)));
  assert(max_diff < 1e-11, sprintf('XtTimesKronSIdTimesY max diff = %e > 1e-11', max_diff));

  results.passed = results.passed + 1;
  fprintf('    [PASS] XtTimesKronSIdTimesY (standard N=100, Ns=4, Nx=25, Ny=25)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('XtTimesKronSIdTimesY standard: %s', err.message);
  fprintf('    [FAIL] XtTimesKronSIdTimesY standard: %s\n', err.message);
end

%% 2. Rectangular / asymmetric dimensions (Nx != Ny)
try
  N = 50;
  Ns = 3;
  Nx = 10;
  Ny = 30;

  x = randn(N*Ns, Nx);
  y = randn(N*Ns, Ny);
  S = randn(Ns, Ns); % asymmetric S

  Z = zeros(N*Ns, Ny);
  for a = 1:Ns
    for b = 1:Ns
      Z((a-1)*N+1:a*N, :) = Z((a-1)*N+1:a*N, :) + S(a,b) * y((b-1)*N+1:b*N, :);
    end
  end
  prod_ref = x' * Z;
  prod_mex = XtTimesKronSIdTimesY(x, S, N, y);

  assert(all(size(prod_mex) == [Nx Ny]), 'Output dimension mismatch');
  max_diff = max(abs(prod_mex(:) - prod_ref(:)));
  assert(max_diff < 1e-11, sprintf('Rectangular test diff = %e > 1e-11', max_diff));

  results.passed = results.passed + 1;
  fprintf('    [PASS] XtTimesKronSIdTimesY (rectangular Nx=10, Ny=30)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('XtTimesKronSIdTimesY rectangular: %s', err.message);
  fprintf('    [FAIL] XtTimesKronSIdTimesY rectangular: %s\n', err.message);
end

%% 3. Edge case: N=1 (single observation per block)
try
  N = 1;
  Ns = 4;
  Nx = 8;
  Ny = 8;

  x = randn(Ns, Nx);
  y = randn(Ns, Ny);
  S = randn(Ns, Ns);

  prod_ref = x' * S * y;
  prod_mex = XtTimesKronSIdTimesY(x, S, N, y);

  max_diff = max(abs(prod_mex(:) - prod_ref(:)));
  assert(max_diff < 1e-12, sprintf('N=1 test diff = %e > 1e-12', max_diff));

  results.passed = results.passed + 1;
  fprintf('    [PASS] XtTimesKronSIdTimesY (edge case N=1)\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('XtTimesKronSIdTimesY N=1: %s', err.message);
  fprintf('    [FAIL] XtTimesKronSIdTimesY N=1: %s\n', err.message);
end

%% 4. Identity covariance S = eye(Ns)
try
  N = 80;
  Ns = 5;
  Nx = 15;
  Ny = 15;

  x = randn(N*Ns, Nx);
  y = randn(N*Ns, Ny);
  S = eye(Ns);

  prod_ref = x' * y;
  prod_mex = XtTimesKronSIdTimesY(x, S, N, y);

  max_diff = max(abs(prod_mex(:) - prod_ref(:)));
  assert(max_diff < 1e-12, sprintf('S=eye test diff = %e > 1e-12', max_diff));

  results.passed = results.passed + 1;
  fprintf('    [PASS] XtTimesKronSIdTimesY (S=eye(Ns))\n');
catch err
  results.failed = results.failed + 1;
  results.errors{end+1} = sprintf('XtTimesKronSIdTimesY S=eye: %s', err.message);
  fprintf('    [FAIL] XtTimesKronSIdTimesY S=eye: %s\n', err.message);
end

end
