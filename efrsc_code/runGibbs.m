clear;
if (isOctave)
  debug_on_error(1);
else
  dbstop if error;
end

global verbosity;
verbosity  = 2; % controls amount of output (0 = some, 1 =
                % more, 2 = most)

%% set initial values -- uses past results
load rsq01.mat; % past results
prefix = 'mainSpec'; % prefix to attach to output files

% compute posterior means
nSaved = numel(chain.theta);
post = floor(nSaved/2):nSaved;
for k=1:numel(chain.beta);
  beta{k} = mean(chain.beta{k}(post,:),1)';
end
%beta{1} = [beta{1}; 0; 0; 0];
%beta{2} = [beta{2}; 0; 0; 0];
Sigma = squeeze(mean(chain.sig(post,:,:),1));
gamma = squeeze(mean(chain.gamma(post,:,:),1));
theta = squeeze(mean(chain.theta(post)));
shape = squeeze(mean(chain.shape(post)));
rho = squeeze(mean(chain.rho(post)));
bll = [mean(chain.bll(post,:),1)'];
sigll = mean(chain.sigll(post));
Sigma = combineSigOPbeta(Sigma(1:2,1:2),gamma,rho,2);
logOmHi = Inf;
varHi = 10;
clear data chain;

%% Load data and set options
make; % compiles *.c into *.mex files if necessary, you will likely need to
      % modify make.m depending on your system
lx{1} = [32:37 40 42 43 44 52 53 54];
lx{2} = lx{1}; lx{3} = lx{1};
lxll = lx{1}; lxhs = lx{1};
data = loadData('al.csv',true,2004,lx,lxll,lxhs);
data.varLo = 0; % lowest sigma_lambda allowed
data.varHi = 4*var(log(1+data.totalSpend(~isnan(data.totalSpend)))); % highest sigma_lambda
data.logOmegaHi = Inf; % highest log omega
data.lamloHi = 3000; % highest kappa


%% create output file with unique name
nameroot=prefix
filename='';
files=dir(sprintf('%s*mat',nameroot));
if (numel(files)==0)
  filename = sprintf('%s01.mat',nameroot);
else 
  mynum = 0;
  for n=1:numel(files)
    j = sscanf(files(n).name,[nameroot '%2d'],1);
    if (j>mynum && files(n).bytes>0)
      mynum = j;
    end
  end
  j = j+1;
  filename = sprintf('%s%02d.mat',nameroot,j);
end
system(sprintf('touch %s',filename));

%% drop select
% sel = squeeze(data.avail(5,2,:)==1);
% data.avail(1,2,sel)= 0;
% data.deduct(1,2,sel)=-999;
% data.maxoop(1,2,sel)=-999;
% data = dropBad(data,data.choice(2,:)'==1 & sel);

% set seed
seed = sum(100*clock);
setSeed(seed);

data.maxIter = 200;
data.verbose = 0;
options.burn = 0; % number of iterations to discard
options.save = 250; % number of iterations to save
options.skip = 50; % save every skipth iteration
options.display = 25; % display info every this many iterations
options.saveInterval = 500; % save a restart file to disk every this many interations
options.savefile = ['tmp_' nameroot '.mat']; % name of restart file
options.hetsked = false; % allow heteroskedasticity?
options.multMH = false; % use multiplicative moral hazard instead of additive?

data.default(isnan(data.choice)) = 0; % NaN's cause problems later
data.choice(isnan(data.choice))=-1;
[status txt] = system('grep processor /proc/cpuinfo | tail -1');
nproc = sscanf(txt(regexp(txt,'\d'):end),'%d')+1;
data.threads=nproc; % number of threads = number of processors
fprintf('Using %d threads.\n',data.threads);
setenv('OMP_NUM_THREADS',num2str(data.threads));

data.prior.beta = [];
data.prior.Sigma=[]; %.A = SigOP;
data.maxIter = 500;
data.verbose = 0;
data.order = 24;

alpha = zeros(size(data.xhs,2),1);

%% drop default people
% data.choice(data.default==1)=-1;
% for i=1:data.N
%   for t=1:data.T
%     if data.default(t,i)==1
%       data.deduct(:,t,i) = NaN;
%       data.totalSpend(t,i) = NaN;
%     end
%   end
% end

data = dropBad(data, ((1:data.N)>200)'); % use small dataset for debugging

%% run mcmc 
[chain data]=gibbs(data, beta, Sigma, gamma, rho, theta, shape, bll,sigll,alpha,options);

% save results
if (isOctave)
  save('-v7',filename,'chain','data','seed');
else
  save(filename,'chain','data','seed');
end
