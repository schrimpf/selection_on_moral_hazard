if exist('mainSpec01.mat','file') || exist('../mainSpec01.mat','file')
  resultFile = 'mainSpec01';
  prefix = 'mainSpec';
else
  resultFile = 'rsq01'; % .mat file of mcmc results
  prefix = 'rsq'; % prefix for tables and figures
end