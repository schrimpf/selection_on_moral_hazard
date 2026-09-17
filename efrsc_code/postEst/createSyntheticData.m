%% createSyntheticData.m
% Loads model results from options/config, generates privacy-preserving synthetic data
% using simulate(..., options.synthetic=true), and saves the synthetic dataset as:
%   1. A .mat file containing the sim struct (with observation fields stripped)
%   2. A .csv file formatted identically to the original al.csv (54 columns with header)
%
% Usage:
%   [sim, csvFile, matFile] = createSyntheticData();
%   [sim, csvFile, matFile] = createSyntheticData(options);
%
% Options fields (all optional):
%   options.resultFile : Results .mat file to load (default: from config.m or revisedCode01)
%   options.nSim / .N  : Number of synthetic individuals (default: data.N)
%   options.csvFile    : Output path for synthetic CSV (default: csv/<prefix>Synthetic.csv)
%   options.matFile    : Output path for synthetic MAT (default: csv/<prefix>Synthetic.mat)
%   options.seed       : RNG seed (default: 1337)
%   options.noiseScale : Continuous noise scaling factor (default: 0.10)
%   options.flipProb   : Categorical randomized response probability (default: 0.05)
%   options.avgPlan    : Whether to average 2004 plan menus (default: true)

function [sim, csvFile, matFile] = createSyntheticData(options)
  if (nargin < 1 || isempty(options))
    options = struct();
  end
  if (~isstruct(options))
    error('createSyntheticData: options must be a struct.');
  end

  % Add parent directory to path to access helper functions and model files
  addpath('..');
  if exist('config.m', 'file')
    config;
  end

  % Determine results file to load
  if isfield(options, 'resultFile')
    resFile = options.resultFile;
  else
    error('no result file specified in options.resultFile');
  end

  if (~exist('prefix', 'var'))
    if isfield(options, 'prefix')
      prefix = options.prefix;
    else
      prefix = 'synthetic';
    end
  end

  fprintf('Loading estimation results from %s...\n', resFile);
  load(resFile);

  % Set random seed
  if isfield(options, 'seed')
    seed = options.seed;
  else
    seed = 1337;
  end
  setSeed(seed);

  % Compute posterior means from chain
  nSaved = numel(chain.meanL);
  if isfield(options, 'post')
    post = options.post;
  else
    post = floor(nSaved/2):nSaved;
  end

  for k = 1:numel(chain.beta)
    beta{k} = mean(chain.beta{k}(post,:), 1)';
  end
  Sigma = squeeze(mean(chain.sig(post,:,:), 1));
  gamma = squeeze(mean(chain.gamma(post,:,:), 1));
  theta = squeeze(mean(chain.theta(post)));
  shape = squeeze(mean(chain.shape(post)));
  rho   = squeeze(mean(chain.rho(post)));
  bll   = mean(chain.bll(post,:), 1)';
  sigll = mean(chain.sigll(post));
  if isfield(chain, 'alpha')
    alpha = mean(chain.alpha(post,:), 1)';
  else
    alpha = 0;
  end

  % Set simulation options
  options.synthetic = true;
  if ~isfield(options, 'avgPlan')
    options.avgPlan = true;
  end

  if isfield(options, 'nSim')
    nSim = options.nSim;
  elseif isfield(options, 'N')
    nSim = options.N;
  else
    nSim = data.N;
  end

  if ~isfield(options, 'sampleAll')
    options.sampleAll = (mod(nSim, data.N) == 0 && nSim >= data.N);
  end

  fprintf('Simulating synthetic dataset with N=%d observations...\n', nSim);
  sim = simulate(nSim, data, beta, Sigma, gamma, rho, shape, theta, bll, sigll, alpha, options);

  % Define output paths
  if ~exist('csv', 'dir')
    mkdir('csv');
  end

  if isfield(options, 'matFile')
    matFile = options.matFile;
  else
    matFile = fullfile('csv', sprintf('%sSynthetic.mat', prefix));
  end

  if isfield(options, 'csvFile')
    csvFile = options.csvFile;
  else
    csvFile = fullfile('csv', sprintf('%sSynthetic.csv', prefix));
  end

  % Save .mat file
  fprintf('Saving synthetic data structure to %s...\n', matFile);
  save(matFile, 'sim');

  % Build 54-column rectangular matrix in exact al.csv format
  fprintf('Exporting synthetic data in al.csv format to %s...\n', csvFile);
  years = unique(sim.year);
  if isempty(years)
    years = [2003; 2004];
  end
  T = numel(years);
  N = sim.N;
  M = zeros(N * T, 54);  % pre-allocate worst case

  row = 0;
  for i = 1:N
    for t = 1:T
      if sim.choice(t, i) <= 0
        continue;  % skip inactive person-years (not at firm this year)
      end
      row = row + 1;
      M(row, 1) = i;                              % 1: id
      M(row, 2) = years(t);                       % 2: year
      M(row, 3) = sim.choice(t, i);               % 3: choice
      M(row, 4:8) = squeeze(sim.avail(1:5, t, i))';   % 4..8: avail1..5
      M(row, 9:13) = squeeze(sim.prem(1:5, t, i))';  % 9..13: price1..5
      M(row, 14:18) = squeeze(sim.deduct(1:5, t, i))'; % 14..18: deduct1..5
      M(row, 19:23) = squeeze(sim.copay(1:5, t, i))';  % 19..23: rate1..5
      M(row, 24:28) = squeeze(sim.maxoop(1:5, t, i))'; % 24..28: max1..5
      M(row, 29) = double(years(t) == 2004);      % 29: d_year_2004
      M(row, 30) = double(years(t) == 2005);      % 30: d_year_2005
      M(row, 31) = double(years(t) == 2006);      % 31: d_year_2006
      M(row, 32:34) = sim.x{3}{t}(i, 2:4);       % 32..34: dcovg_tierr_2..4
      M(row, 35:37) = sim.x{3}{t}(i, 5:7);       % 35..37: zgroup_2..4
      M(row, 38) = 1 - sum(sim.x{3}{t}(i, 5:7)); % 38: zgroup_5
      M(row, 39) = sim.totalSpend(t, i);          % 39: total
      age_it = sim.x{3}{t}(i, 8);
      M(row, 40) = age_it;                       % 40: age
      M(row, 41) = (age_it^2) / 100;              % 41: age2
      M(row, 42) = sim.x{3}{t}(i, 9);            % 42: sex
      M(row, 43) = sim.x{3}{t}(i, 10);           % 43: tenure_hr
      wage_it = sim.x{3}{t}(i, 11);
      M(row, 44) = wage_it;                      % 44: avgWage
      M(row, 45) = (wage_it^2) / 100;             % 45: wage2
      if (sim.x{3}{t}(i, 2) == 1)
        M(row, 46) = 2;                          % 46: family_size
      elseif (sim.x{3}{t}(i, 4) == 1)
        M(row, 46) = 1;
      else
        M(row, 46) = 0;
      end
      M(row, 47) = 1.0;                           % 47: rs_fam
      M(row, 48) = 1.0;                           % 48: rs_emp
      M(row, 49) = 0.0;                           % 49: logrsfam
      M(row, 50) = 0;                             % 50: default
      M(row, 51) = 1 - sum(sim.x{3}{t}(i, 12:14)); % 51: rsfq1
      M(row, 52:54) = sim.x{3}{t}(i, 12:14);     % 52..54: rsfq2..4
    end
  end
  M = M(1:row, :);  % trim to active rows only

  header = 'id,year,choice,avail1,avail2,avail3,avail4,avail5,price1,price2,price3,price4,price5,deduct1,deduct2,deduct3,deduct4,deduct5,rate1,rate2,rate3,rate4,rate5,max1,max2,max3,max4,max5,d_year_2004,d_year_2005,d_year_2006,dcovg_tierr_2,dcovg_tierr_3,dcovg_tierr_4,zgroup_2,zgroup_3,zgroup_4,zgroup_5,total,age,age2,sex,tenure_hr,avgWage,wage2,family_size,rs_fam,rs_emp,logrsfam,default,rsfq1,rsfq2,rsfq3,rsfq4';

  fid = fopen(csvFile, 'w');
  if fid < 0
    error('Could not open %s for writing.', csvFile);
  end
  fprintf(fid, '%s\n', header);
  fclose(fid);
  dlmwrite(csvFile, M, '-append', 'delimiter', ',');

  fprintf('Successfully created synthetic data:\n  MAT: %s\n  CSV: %s\n', matFile, csvFile);
end
