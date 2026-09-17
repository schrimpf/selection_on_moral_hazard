function success = run_test_suite()
%% run_test_suite.m
% Master test runner for the insurance choice Gibbs estimation pipeline.
% Discovers and executes all test suites, reporting aggregated results.

close all;
clear all;

% Ensure path includes project root and tests directory
testsDir = fileparts(mfilename('fullpath'));
if isempty(testsDir)
  testsDir = pwd;
end
addpath(testsDir);
addpath(fullfile(testsDir, '..'));

global verbosity;
verbosity = 0;

suites = {
  'test_math_utils',       @test_math_utils;
  'test_matrix_kernels',   @test_matrix_kernels;
  'test_choice_engine',    @test_choice_engine;
  'test_latent_samplers',  @test_latent_samplers;
  'test_gibbs_pipeline',   @test_gibbs_pipeline;
  'test_sample_mu',        @test_sample_mu;
};

fprintf('\n======================================================================\n');
fprintf('RUNNING GIBBS ESTIMATION TEST SUITE\n');
fprintf('======================================================================\n\n');

total_passed = 0;
total_failed = 0;
all_errors = {};
suite_results = struct();

overall_tic = tic;

for s = 1:size(suites, 1)
  suite_name = suites{s, 1};
  suite_func = suites{s, 2};

  t0 = tic;
  try
    res = suite_func();
    duration = toc(t0);
    suite_results.(suite_name).status = (res.failed == 0);
    suite_results.(suite_name).duration = duration;
    suite_results.(suite_name).passed = res.passed;
    suite_results.(suite_name).failed = res.failed;

    total_passed = total_passed + res.passed;
    total_failed = total_failed + res.failed;
    if ~isempty(res.errors)
      all_errors = [all_errors, res.errors];
    end
  catch err
    duration = toc(t0);
    suite_results.(suite_name).status = false;
    suite_results.(suite_name).duration = duration;
    suite_results.(suite_name).passed = 0;
    suite_results.(suite_name).failed = 1;
    total_failed = total_failed + 1;
    all_errors{end+1} = sprintf('%s suite crash: %s', suite_name, err.message);
    fprintf('    [CRASH] %s: %s\n', suite_name, err.message);
  end
  fprintf('\n');
end

overall_duration = toc(overall_tic);

%% Print Summary Table
fprintf('======================================================================\n');
fprintf('TEST SUITE SUMMARY\n');
fprintf('======================================================================\n');
fprintf('%-25s | %8s | %8s | %8s | %10s\n', 'Test Suite', 'Passed', 'Failed', 'Status', 'Time (s)');
fprintf('----------------------------------------------------------------------\n');

all_suites_passed = true;
for s = 1:size(suites, 1)
  name = suites{s, 1};
  info = suite_results.(name);
  if info.status
    status_str = 'PASSED';
  else
    status_str = 'FAILED';
    all_suites_passed = false;
  end
  fprintf('%-25s | %8d | %8d | %8s | %10.3f\n', ...
          name, info.passed, info.failed, status_str, info.duration);
end

fprintf('======================================================================\n');
if all_suites_passed
  fprintf('OVERALL RESULT: SUCCESS (%d checks passed in %.2f s)\n', total_passed, overall_duration);
else
  fprintf('OVERALL RESULT: FAILURE (%d passed, %d failed in %.2f s)\n', ...
          total_passed, total_failed, overall_duration);
  fprintf('\nFailed checks:\n');
  for e = 1:numel(all_errors)
    fprintf('  - %s\n', all_errors{e});
  end
end
fprintf('======================================================================\n\n');

success = all_suites_passed;
end
