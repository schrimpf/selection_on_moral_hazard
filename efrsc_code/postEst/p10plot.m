close all;
path(path,'..');
config;
load(resultFile);
nSim = data.N*10;
seed = 1337; %sum(100*clock);
setSeed(seed);
clear options;
clear b15data;
options.cf = true;
options.nosel = true;
options.sampleAll=true;
options.avgPlan = true;
options.multMH = false;

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
if isfield(chain,'alpha')
  alpha = mean(chain.alpha(post,:),1)';
else
  alpha = 0;
end


names = {'\omega','\psi','\mu_{\lambda 1}','\mu_{\lambda 2}', ...
         '\mu_{\lambda 3}'};
xnames{1} = {'Constant','family','+ spouse','+ child'}'
xnames{2} = xnames{1};
xnames{3} = {'Constant','family','+ spouse','+ child', ...
             'switch 04','switch 05','switch 06','2004','2005'};
[order xint wint] = gqzero(data.order);
xint = xint'*sqrt(2);
order = numel(xint);
wint = wint/sqrt(pi);

%% 5-1 90-10 data
setSeed(seed);
options.balance15 = 0.5;
b15data{1}=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                    options);
% spending changes
dsFN = b15data{1}.spendFI(2,:) - b15data{1}.spendNI(2,:);
s5 = b15data{1};
s5.choice(:) = 5;
s5.totalSpend = spending(s5);
s1 = b15data{1};
s1.choice(:) = 1;
s1.totalSpend = spending(s1);
ds = s5.totalSpend(2,:) - s1.totalSpend(2,:);


options.avgPlan = true;
portion = 0.05:0.05:0.95
for i=1:numel(portion)
  setSeed(seed);
  options.balance15 = portion(i);
  b15data{i}=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                      options);
  epsi(i) = mean(log(b15data{i}.psi(b15data{i}.choice(2,:)==1)));
  eomega(i) = mean(log(b15data{i}.omega(b15data{i}.choice(2,:)==1)));
  emul(i) = mean(b15data{i}.muL(b15data{i}.choice(2,:)==1,2));
  eds(i) = mean(ds(b15data{i}.choice(2,:)==1 & ~isnan(ds)));
end

figure;
subplot(2,2,1),plot(portion,eds)
title('E[\Delta spend 5-1|choose 1]');
xlabel('% choose 5');
subplot(2,2,2),plot(portion,eomega)
title('E[log \omega|choose 1]');
xlabel('% choose 5');
subplot(2,2,3),plot(portion,emul)
title('E[\lambda|choose 1]');
xlabel('% choose 5');
subplot(2,2,4),plot(portion,epsi);
title('E[log \psi|choose 1]');
xlabel('% choose 5');
print('-depsc2',sprintf('figures/%sEcond1',prefix));

figure;
subplot(2,2,1),plot(1-portion,eds)
title('E[\Delta spend 5-1|choose 1]');
xlabel('% choose 1');
subplot(2,2,2),plot(1-portion,eomega)
title('E[log \omega|choose 1]');
xlabel('% choose 1');
subplot(2,2,3),plot(1-portion,emul)
title('E[\mu_\lambda|choose 1]');
xlabel('% choose 1');
subplot(2,2,4),plot(1-portion,epsi);
title('E[log \psi|choose 1]');
xlabel('% choose 1');
print('-depsc2',sprintf('figures/%sEcond1v2',prefix));
writecsv(sprintf('csv/%sselectMH.csv',prefix), ...
         {'portion 1','change spend 5-1', ...
          'log omega','mu_lambda','log psi'}, ...
         [1-portion', eds', eomega', emul', epsi']);
