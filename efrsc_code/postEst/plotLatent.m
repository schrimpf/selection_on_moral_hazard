clear;
path('../ver2',path);
path('../ver2/counterfactuals',path);
config;
load(resultFile);
seed = 1337;
setSeed(seed);
nSim = data.N*10;

%% compute posterior means
nSaved = numel(chain.meanL);
post = floor(nSaved/2):nSaved;
for k=1:numel(chain.beta);
  beta{k} = mean(chain.beta{k}(post,:),1)';
end
Sigma = squeeze(mean(chain.sig(post,:,:),1));
gamma = squeeze(mean(chain.gamma(post,:,:),1));
theta = squeeze(mean(chain.theta(post)));
shape = squeeze(mean(chain.shape(post)));
rho = squeeze(mean(chain.rho(post)));
bll = mean(chain.bll(post,:),1)';
sigll = mean(chain.sigll(post));

%% plot risk scores and spending
sobs = data.totalSpend(2,:)';
figure;
subplot(2,1,1),plot(log(data.x{3}{2}(~isnan(sobs),end-1)),log(1+sobs(~isnan(sobs))),'.')
xlabel('log risk score');
ylabel('log(1+spend)');
subplot(2,1,2),plot((data.x{3}{2}(~isnan(sobs),end-1)),log(1+sobs(~isnan(sobs))),'.')
xlabel('risk score');
ylabel('log(1+spend)');
print('-dpng','rsSpend.png');


%% show that high risk scores cause spending outliers
options.sampleAll = true;
tic
simdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll, ...
                 options);
toc
%% plot risk scores and spending
ssim = simdata.totalSpend(2,:)';
figure;
subplot(2,1,1),plot(log(simdata.x{3}{2}(~isnan(ssim),end-1)), ...
                    log(1+ssim(~isnan(ssim))),'c.', ...
                    log(data.x{3}{2}(~isnan(sobs),end-1)),log(1+sobs(~isnan(sobs))),'r.');
xlabel('log risk score');
ylabel('log(1+spend)');
legend('simulated','observed');
subplot(2,1,2),plot((simdata.x{3}{2}(~isnan(ssim),end-1)), ...
                    log(1+ssim(~isnan(ssim))),'c.', ...
                    (data.x{3}{2}(~isnan(sobs),end-1)),log(1+sobs(~isnan(sobs))),'r.');
xlabel('risk score');
ylabel('log(1+spend)');
legend('simulated','observed');
print('-dpng','rsSpendSim.png');

rs = data.x{3}{2}(:,end-1);
rss = simdata.x{3}{2}(:,end-1);
[mean(ssim(~isnan(ssim))) mean(sobs(~isnan(sobs)))]
[mean(ssim(~isnan(ssim) & rss<4)) mean(sobs(~isnan(sobs) & rs<4))]
[mean(ssim(~isnan(ssim) & rss>4)) mean(sobs(~isnan(sobs) & rs>4))]

