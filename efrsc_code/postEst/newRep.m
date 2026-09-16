clear; close all;
config;
path(path,'..');
load(resultFile);
negLambdaZeroUtil = false;
nSim = data.N*10;
seed = 1337; %sum(100*clock);
setSeed(seed);

options.cf = true;
options.nosel = false;
options.sampleAll=true;
options.multMH = false;
% compute posterior means
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
xnames{1} = {'Constant','family','+ spouse','+ child', ...
             'switch 04','switch 05','switch 06','switch later'};
xnames{2} = xnames{1};
xnames{3} = {'Constant','family','+ spouse','+ child', ...
             'switch 04','switch 05','switch 06','switch later','2004','2005'};
[order xint wint] = gqzero(data.order);
xint = xint'*sqrt(2);
order = numel(xint);
wint = wint/sqrt(pi);

%%
simdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha,options);

%% spending stats and densities
nbin = 64;
path(path,'kde');
maxs = max(log(data.totalSpend(:)+1));
select = squeeze(data.avail(4,:,:)==1);
subset =  ~select(:) & ~isnan(data.totalSpend(:));
simsub = ~squeeze(simdata.avail(4,:,:)==1) & ~isnan(simdata.totalSpend);
simsub = simsub(:);
out = fopen(sprintf('tex/tables/%sfit.tex',prefix),'w');
fprintf(out,['& Observed & Estimated \\\\ \\hline \n' ...
             '\\multicolumn{3}{c}{Flex} \\\\ \\hline \n']);
for(c=1:3) 
  fprintf(out,'$P(c=%d)$ & %.2f & %.2f \\\\ \n',c, mean(data.choice(subset)==c), ...
          mean(simdata.choice(simsub)==c));
end
fprintf(out,'Mean spending & %.0f & %.0f \\\\ \n', mean(data.totalSpend(subset)), ...
         mean(simdata.totalSpend(simsub)));
fprintf(out,'Median spending & %.0f & %.0f \\\\ \n', median(data.totalSpend(subset)), ...
         median(simdata.totalSpend(simsub)));
fprintf(out,'$P($spend$=0)$ & %.2f & %.2f \\\\ \n', mean(data.totalSpend(subset)==0), ...
         mean(simdata.totalSpend(simsub)==0));

figure; colormap('gray');
%[bw dens x] = kde(log(data.totalSpend(subset)+1),nbin,0,maxs);
%[bw fdens xf] = kde(log(simdata.totalSpend(simsub)+1),nbin,0,maxs);
nbin = 20;
bin = [-inf 0.5:0.5:11 inf];
[dens] = histc(log(data.totalSpend(subset)+1),bin);
[fdens] = histc(log(simdata.totalSpend(simsub)+1),bin);
x = bin';
x(1) = 0;
x(end) = 11.5;
xf = x;
dens = dens/sum(dens);
fdens = fdens/sum(fdens);
%plot(x,dens,'r-',xf,fdens,'b--');
bar(x,[dens fdens]);
%text(.1,0.2,'P(s=0)=%.2f',mean(data.totalSpend(subset)==0),'Color','r');
%text(.1,0.1,'P(s=0)=%.2f',mean(simdata.totalSpend(simsub)==0),'Color','b');
%ylim([0 0.35]);
legend('Observed','Estimated');
xlabel('log dollars');
ylabel('density');
print('-depsc2',sprintf('figures/%sSpendFlex',prefix));
writecsv(sprintf('csv/%sspendFlex.csv',prefix), ...
         {'log spend','observed','fitted'}, ...
         [x',dens',fdens']);



select = squeeze(data.avail(4,:,:)==1);
subset =  select(:) & ~isnan(data.totalSpend(:));
simsub = squeeze(simdata.avail(4,:,:)==1) & ~isnan(simdata.totalSpend);
simsub = simsub(:);
fprintf(out,'\\multicolumn{3}{c}{Select} \\\\ \\hline \n');
for(c=1:5) 
  fprintf(out,'$P(c=%d)$ & %.2f & %.2f \\\\ \n',c, mean(data.choice(subset)==c), ...
          mean(simdata.choice(simsub)==c));
end
fprintf(out,'Mean spending & %.0f & %.0f \\\\ \n', mean(data.totalSpend(subset)), ...
        mean(simdata.totalSpend(simsub)));
fprintf(out,'Median spending & %.0f & %.0f \\\\ \n', median(data.totalSpend(subset)), ...
        median(simdata.totalSpend(simsub)));
fprintf(out,'$P($spend$=0)$ & %.2f & %.2f \\\\ \n', mean(data.totalSpend(subset)==0), ...
        mean(simdata.totalSpend(simsub)==0));
fprintf(out,'\\hline ');

figure; colormap('gray');
%[bw dens x] = kde(log(data.totalSpend(subset)+1),nbin,0,maxs);
%[bw fdens xf] = kde(log(simdata.totalSpend(simsub)+1),nbin,0,maxs);
[dens] = histc(log(data.totalSpend(subset)+1),bin);
[fdens] = histc(log(simdata.totalSpend(simsub)+1),bin);
dens = dens/sum(dens);
fdens = fdens/sum(fdens);
bar(x,[dens fdens]);
%plot(x,dens,'r-',xf,fdens,'b--');
%ylim([0 0.35]);
%text(.1,0.2,'P(s=0)=%.2f',mean(data.totalSpend(subset)==0),'Color','r');
%text(.1,0.1,'P(s=0)=%.2f',mean(simdata.totalSpend(simsub)==0),'Color','b');
legend('Observed','Estimated');
xlabel('log dollars');
ylabel('density');
print('-depsc2',sprintf('figures/%sSpendSel',prefix));
writecsv(sprintf('csv/%sspendSelect.csv',prefix), ...
         {'log spend','observed','fitted'}, ...
         [x',dens',fdens']);

nonsingle = sum(data.x{1}(:,2:4),2);
subset = ~isnan(data.totalSpend) & ~(ones(data.T,1)*nonsingle');
nonsingleSim = sum(simdata.x{1}(:,2:4),2);
simsub = ~isnan(simdata.totalSpend) & ~(ones(data.T,1)*nonsingleSim');
fprintf(out,'\\hline \\multicolumn{3}{c}{Single (N = %d)} \\\\ \n',sum(subset(:)));
fprintf(out,'Mean spending & %.0f & %.0f \\\\ \n', mean(data.totalSpend(subset)), ...
         mean(simdata.totalSpend(simsub)));
fprintf(out,'Median spending & %.0f & %.0f \\\\ \n', median(data.totalSpend(subset)), ...
         median(simdata.totalSpend(simsub)));
fprintf(out,'P(spending=0) & %.2f & %.2f \\\\ \n', mean(data.totalSpend(subset)==0), ...
         mean(simdata.totalSpend(simsub)==0));
subset = ~isnan(data.totalSpend) & ones(data.T,1)*nonsingle';
simsub = ~isnan(simdata.totalSpend) & ones(data.T,1)*nonsingleSim';
fprintf(out,'\\hline \\multicolumn{3}{c}{Non-single (N = %d)} \\\\ \n',sum(subset(:)));
fprintf(out,'Mean spending & %.0f & %.0f \\\\ \n', mean(data.totalSpend(subset)), ...
         mean(simdata.totalSpend(simsub)));
fprintf(out,'Median spending & %.0f & %.0f \\\\ \n', median(data.totalSpend(subset)), ...
         median(simdata.totalSpend(simsub)));
fprintf(out,'P(spending=0) & %.2f & %.2f \\\\ \n', mean(data.totalSpend(subset)==0), ...
         mean(simdata.totalSpend(simsub)==0));
fclose(out);

%% effect of select plans
%options.sampleCondChoice = true;
%cdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll, ...
%               options);
options.sampleCondChoice = false;
options.avgPlan = true;
avgdata=simulate(nSim,data,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
out = fopen(sprintf('tex/tables/%sSpendEff.tex',prefix),'w');
fprintf(out,[' & \\multicolumn{2}{c}{Flex-Select} & \\multicolumn{2}{c}{Full-No} ' ...
             '\\\\ \n' ...
             'Sample & \\$ & log \\$ & \\$ & log \\$ \\\\ \\hline \n']);
%dat = {simdata, cdata, avgdata};
dat = {simdata, avgdata};
%name = {'Select 04','Select 04 | obs choice','All 04'};
name = {'Select 04','All 04'};
for d=1:numel(dat)
  ds = dat{d}.spendNoSel - dat{d}.totalSpend ;
  select = squeeze(dat{d}.avail(4,2,:)==1);
  ind = ~isnan(ds(2,:)) & select';
  dsl = log(1+dat{d}.spendNoSel) - log(1+dat{d}.totalSpend);
  dsfn = dat{d}.spendFI - dat{d}.spendNI;
  %ifn = ~isnan(dsfn(2,:));
  dslfn = log(dat{d}.spendFI+1) - log(dat{d}.spendNI+1);
  fprintf(out,'%s & %.1f & %.3g & %.1f & %.3g \\\\ \n', name{d}, ...
          mean(ds(2,ind)), mean(dsl(2,ind)), mean(dsfn(2,ind)), ...
          mean(dslfn(2,ind)));
end
fclose(out);

%% 5-1 90-10 data
setSeed(seed);

% load data without plan 1 deleted
if (negLambdaZeroUtil)
  data2 = loadData('al.csv',true,2004);
  data2.order = 12;
  data2.varLo = 0;
  data2.varHi = 4*var(log(1+data2.totalSpend(~isnan(data2.totalSpend))));
  avgSpend = mean(data2.totalSpend(~isnan(data2.totalSpend)));
  data2.logOmegaHi = Inf;
  data2.lamloHi = 3000;
  data2.maxIter = data.maxIter;
  data2.order = data.order;
  data2.threads = data.threads;
  sel = squeeze(data2.avail(5,2,:)==1);
  data2 = dropBad(data2,data2.choice(2,:)'==1 & sel);
else 
  data2 = data;
end
nSim = 10*data2.N;

options.avgPlan = true;
options.balance15 = 0.9;
b15data=simulate(nSim,data2,beta,Sigma,gamma,rho,shape,theta,bll,sigll, alpha,...
                 options);
%% density of spending changes 
ds = b15data.spendFI(2,:) - b15data.spendNI(2,:);
figure;
hi = 2000;
[bw dens x] = kde(ds(~isnan(ds)),nbin,0,hi);
dsFN = ds;
plot(x,dens);
xlabel('dollars');
ylabel('density');
print('-depsc2',sprintf('figures/%sdsFN',prefix));
figure;
s5 = b15data;
s5.choice(:) = 5;
s5.totalSpend = spending(s5);
s1 = b15data;
s1.choice(:) = 1;
s1.totalSpend = spending(s1);
ds = s5.totalSpend(2,:) - s1.totalSpend(2,:);
[bw dens x] = kde(ds(~isnan(ds)),nbin,0,hi);
plot(x,dens);
xlabel('dollars');
ylabel('density');
print('-depsc2',sprintf('figures/%sds51',prefix));


options.balance15 = 0.5;
tmpdata1=simulate(nSim,data2,beta,Sigma,gamma,rho,shape,theta,bll,sigll, alpha, ...
                 options);
setSeed(seed);

options.balance15 = 0.1;
tmpdata2=simulate(nSim,data2,beta,Sigma,gamma,rho,shape,theta,bll,sigll,alpha, ...
                 options);
options.balance15 = 0.9;
out = fopen(sprintf('tex/tables/%sSpendP5.tex',prefix),'w');
fprintf(out,'& 0\\%% & 10\\%% & 50\\%% & 90\\%% & 100\\%% \\\\ \n $E[spend]$');
fn = @(x) mean(x.totalSpend(2,~isnan(x.totalSpend(2,:))));
fprintf(out,'&  %.0f ',[fn(s1) fn(tmpdata2) fn(tmpdata1) fn(b15data) fn(s5)]);
fprintf(out,' \\\\ \n $E[\\Delta_{5-1} spend | choose\\, 1]$ ');
fn = @(x) mean(ds(x.choice(2,:)==1));
fprintf(out, '& %.0f ',[fn(s1),fn(tmpdata2),fn(tmpdata1),fn(b15data)]);
fprintf(out,' & -- \\\\ \n');
fprintf(out,' \\\\ \n $E[\\Delta_{F-N} spend | choose\\, 1]$ ');
fn = @(x) mean(dsFN(x.choice(2,:)==1));
fprintf(out, '& %.0f ',[fn(s1),fn(tmpdata2),fn(tmpdata1),fn(b15data)]);
fprintf(out,' & -- \\\\ \n');
fprintf(out,' \\\\ \n $Med[\\Delta_{5-1} spend | choose\\, 1]$ ');
fn = @(x) median(ds(x.choice(2,:)==1));
fprintf(out, '& %.0f ',[fn(s1),fn(tmpdata2),fn(tmpdata1),fn(b15data)]);
fprintf(out,' & -- \\\\ \n');
fprintf(out,' \\\\ \n $Med[\\Delta_{F-N} spend | choose\\, 1]$ ');
fn = @(x) median(dsFN(x.choice(2,:)==1));
fprintf(out, '& %.0f ',[fn(s1),fn(tmpdata2),fn(tmpdata1),fn(b15data)]);
fprintf(out,' & -- \\\\ \n');
fclose(out);

out = fopen(sprintf('tex/tables/%sds.tex',prefix),'w');
tau = [0.1 0.25 0.50 0.75 0.9];
fprintf(out,' & Mean & Std Dev ');
fprintf(out,'& %.2f ',tau);
fprintf(out,' & Mean$\\vert$choose 1 \\\\ \\hline \n');
y = ds(~isnan(ds));
fprintf(out,' 5-1 '); 
fprintf(out,' & %.0f ',[mean(y) std(y) quantile(y',tau)' mean(ds(b15data.choice(2,:)==1))]);
fprintf(out,' \\\\ \n Full-No ');
y = dsFN(~isnan(dsFN));
fprintf(out,' & %.0f ',[mean(y) std(y) quantile(y',tau)' mean(dsFN(b15data.choice(2,:)==1))]);
fprintf(out, ' \\\\ \n');
fclose(out);

% spending effects and omega by x's
TElev = (simdata.spendNoSel(2,:) - simdata.totalSpend(2,:))' ;
select = squeeze(simdata.avail(4,2,:)==1);
ind = ~isnan(TElev) & select;
TElog = (log(1+simdata.spendNoSel(2,:)) - ...
         log(1+simdata.totalSpend(2,:)))';
out = fopen(sprintf('tex/tables/%sMHbyX.tex',prefix),'w');
fprintf(out,['& & Obs. & Mean $\\omega$ & Mean spend 5-1 & Mean TE ' ...
             '(levels) & Mean TE (logs) \\\\ & \n']);
fprintf(out,' & (%d) ', 1:5); fprintf(out,'\\\\ \n');
sub = simdata.x{3}{2}(:,8)>=43;
fprintf(out,['\\multirow{2}{*}{(1)} & Above median age ($>=$43) &' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
        sum(data.x{3}{2}(:,8)>=43), mean(simdata.omega(sub)), ...
        mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
        mean(TElog(sub & ind)) );       
sub = simdata.x{3}{2}(:,8)<43;
fprintf(out,[' & Below median age ($<$43) &' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
        sum(data.x{3}{2}(:,8)<43), mean(simdata.omega(sub)), ...
        mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
        mean(TElog(sub & ind)) );
sub = simdata.x{3}{2}(:,9)==0;
fprintf(out,['\\multirow{2}{*}{(2)} & Male & ', ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.x{3}{2}(:,9)==0), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
sub = simdata.x{3}{2}(:,9)==1;             
fprintf(out,[' & Female & ', ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.x{3}{2}(:,9)==1), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
sub = simdata.x{3}{2}(:,11) >= 36;
fprintf(out,['\\multirow{2}{*}{(3)} & Above median income ($>=$\\$37,000) & ' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.x{3}{2}(:,11)>=37), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
sub = simdata.x{3}{2}(:,11) < 37;
fprintf(out,[' & Below median income ($<$\\$37,000) & ' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.x{3}{2}(:,11)<37), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
sub = (simdata.choice(1,:)==1 | simdata.choice(1,:)==2)';
fprintf(out,['\\multirow{2}{*}{(4)} & Less coverage in 2003 & ' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.choice(1,:)==1 | data.choice(1,:)==2), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
sub = (simdata.choice(1,:)==3)';
fprintf(out,[' & More coverage in 2003 & ' ...
             ' %d & %.3g & %.3g & %.3g & %.3g  \\\\ \n'], ...
             sum(data.choice(1,:)==3), mean(simdata.omega(sub)), ...
             mean(ds(sub & ~isnan(ds)')), mean(TElev(sub & ind)), ...
             mean(TElog(sub & ind)) );
fclose(out);

%% plots of P(5)
%patho(path,'counterfactuals');
slideTables2(b15data,prefix,sprintf('csv/%sP5latent.csv',prefix));
plotP5psiOm(b15data,prefix);

%% with average x's, no corr
options.avgX = true;
setSeed(seed);

S = Sigma;
S(3:4,1:2) = 0; S(1:2,3:4) = 0; S(2,[1 3:4]) = 0; S([1 3:4],2) = 0; 
S(1,2:4) = 0; S(2:4,1) = 0;
nocdata=simulate(nSim,data2,beta,S,gamma,rho,shape,theta,bll,sigll,alpha,options);
slideTables2(nocdata,[prefix 'nocAX'],sprintf('csv/%sP5latentAvgXnoCorr.csv',prefix));
plotP5psiOm(nocdata,[prefix 'nocAX']);

%% with uncorrelated x's
options.uncorrX = true;
options.avgX = false;
setSeed(seed);

S = Sigma;
S(3:4,1:2) = 0; S(1:2,3:4) = 0; S(2,[1 3:4]) = 0; S([1 3:4],2) = 0; 
S(1,2:4) = 0; S(2:4,1) = 0;
nocdata=simulate(nSim,data2,beta,S,gamma,rho,shape,theta,bll,sigll,alpha,options);
slideTables2(nocdata,[prefix 'uncorrX'],sprintf('csv/%sP5latentuncorr.csv',prefix));
plotP5psiOm(nocdata,[prefix 'uncorrX']);
