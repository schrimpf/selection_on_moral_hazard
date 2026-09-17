config;
if ~exist('force', 'var')
  force = false;
end
path(path,'..');

%% create output folders if needed
folders = {'csv','tex','figures','tex/tables'};
for i=1:numel(folders)
  if exist(folders{i})~=7
    mkdir(folders{i});
  end
end

% helper function to check if target needs rebuild
function need = needsRun(targetFile, resFile, force)
  if force, need = true; return; end
  targetInfo = dir(targetFile);
  if isempty(targetInfo), need = true; return; end
  rInfo = dir(resFile);
  if isempty(rInfo)
    rInfo = dir(['../' resFile]);
  end
  if isempty(rInfo)
    need = false;
  else
    need = (targetInfo(1).datenum < rInfo(1).datenum);
  end
end

if needsRun(sprintf('tex/tables/%sParmP.tex',prefix), [resultFile '.mat'], force)
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING PAPERTABLES\n']);
  paperTables;
  paperTablesCSV;
end

if needsRun(sprintf('csv/%sselectMH.csv',prefix), [resultFile '.mat'], force)
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING P10PLOT\n']);
  p10plot;
end

%fitReport;
if needsRun(sprintf('tex/tables/%sds.csv',prefix), [resultFile '.mat'], force)
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING NEWREP\n']);
  newRep;
end

if needsRun(sprintf('tex/tables/%sWelfareP.csv',prefix), [resultFile '.mat'], force)
  fprintf(['\n############################################################' ...
           '####################\n' ...
           'RUNING WELFARE\n']);
  welfare;
end
