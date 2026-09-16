config;
force = false;
path(path,'..');

%% create output folders if needed
folders = {'csv','tex','figures','tex/tables'};
for i=1:numel(folders)
  if exist(folders{i})~=7
    mkdir(folders{i});
  end
end

% run counterfactuals that haven't been run since resultFile last changed
resInfo = dir([resultFile '.mat']);
check = dir(sprintf('tex/tables/%sParmP.tex',prefix));
if isempty(check) || check.datenum<resInfo.datenum || force
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING PAPERTABLES\n']);
  paperTables;
  paperTablesCSV;
end

resInfo = dir([resultFile '.mat']);
check = dir(sprintf('csv/%sselectMH.csv',prefix));
if isempty(check) || check.datenum<resInfo.datenum || force
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING P10PLOT\n']);
  p10plot;
end


%fitReport;
resInfo = dir([resultFile '.mat']);
check = dir(sprintf('tex/tables/%sds.csv',prefix));
if isempty(check) || check.datenum<resInfo.datenum || force
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING NEWREP\n']);
  newRep;
end
resInfo = dir([resultFile '.mat']);
check = dir(sprintf('tex/tables/%sWelfareP.csv',prefix));
if isempty(check) || check.datenum<resInfo.datenum || force
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING WELFARE\n']);
  welfare;
end
