%% Delete observations from data
% data = data structure produced by loadData.m
% bad = data.N x 1 logical array. Observations with bad==1 will be deleted
function [data bad]=dropBad(data,bad)
  if (sum(bad)==0) 
    return;
  end
  
  for field = {'choice','totalSpend','default'}
    if isfield(data,field{1})
      data.(field{1}) = data.(field{1})(:,~bad);
    end
  end
  for field= {'avail','prem','deduct','copay','maxoop'}
    data.(field{1}) = data.(field{1})(:,:,~bad);
  end
  for j=1:2
    data.x{j} = data.x{j}(~bad,:);
  end
  for t=1:data.T 
    data.x{3}{t} = data.x{3}{t}(~bad,:);
  end
  data.xll = data.xll(~bad,:);
  if isfield(data,'xhs')
    data.xhs = data.xhs(~bad,:);
  end
  data.N = sum(~bad);
  if (isfield(data,'dropped'))
    data.dropped = [data.dropped; data.uid(find(bad))'];
  else 
    data.dropped = data.uid(find(bad))';
    if (size(data.dropped,2)>1) 
      data.dropped = data.dropped';
    end
  end
  data.uid = data.uid(~bad);
  k = [size(data.x{1},2) size(data.x{2},2) size(data.x{3}{1},2)];
  n = size(data.x{1},1);
  data.X = [ data.x{1}    zeros(n,k(2)) zeros(n,k(3)); ...
             zeros(n,k(1))   data.x{2}   zeros(n,k(3));];
  for t=1:data.T
    data.X = [data.X; zeros(n,k(1)+k(2)) data.x{3}{t}];
  end

end