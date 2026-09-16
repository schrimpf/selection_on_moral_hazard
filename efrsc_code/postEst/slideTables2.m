function slideTables2(sdata,prefix,csv)
  if (nargin<3)
    csv='';
  end
%% tables and figures about moral hazard
  subset{1} = @(d) d.choice(:)>0;
  subset{2} = @(d) d.choice(:)>0 & d.select(:)~=1;
  subset{3} = @(d) d.choice(:)>0 & d.select(:)==1;
  subnames = {'All','Flex choice set','Select choice set'};
  nbin = 64;
  
  %% graphs of spending diff and P(5) as functions of omega, etc
  nq = 25;
  bw = 0.05
  ind = sdata.choice(2,:)>0;
  y = [sdata.choice(2,ind)==5]';
  x = squeeze([sdata.exSpend(5,2,ind)-sdata.exSpend(1,2,ind)]);
  figure;
  [q fy]=plotSpend(x,y,nq,@(x) mean(x),true,true,bw);
  xlabel('E[spend|c=5]-E[spend|c=1]');
  ylabel('P(c=5)');
  print('-depsc2',sprintf('figures/%sP5ES51',prefix)); 

  x = sdata.eMoralHaz(ind);
  figure;
  [q fy]=plotSpend(x,y,nq,@(x) mean(x),true,true);
  xlabel('E[spend|full ins]-E[spend|no ins]');
  ylabel('P(c=5)');
  print('-depsc2',sprintf('figures/%sP5ESFI',prefix)); 

  figure;
  y = (sdata.choice(2,ind)==5)'*ones(1,2);
  x = [squeeze([sdata.exSpend(5,2,ind)-sdata.exSpend(1,2,ind)]) ...
       sdata.eMoralHaz(ind)];
  [q fy]=plotSpend(x,y,nq,@(x) mean(x),false,true);
  legend('E[spend|c=5]-E[spend|c=1]','E[spend|full ins]-E[spend|no ins]',0);
  ylabel('P(c=5)');
  xlabel('Quantiles');
  print('-depsc2',sprintf('figures/%sP5ESquant',prefix)); 
  
  y = (sdata.choice(2,ind)==5)'*ones(1,3);
  x = [sdata.omega(ind) sdata.psi(ind) sdata.eLambda(2,ind)'];
  figure;
  [q fy]=plotSpend(x,y,nq,@(x) mean(x),false,true);
  xlabel('Quantiles');
  ylabel('P(c=5)');
  legend('\omega','\psi','E[\lambda]',0);
  print('-depsc2',sprintf('figures/%sP5latent',prefix)); 
  out = fopen(sprintf('tex/tables/%sP5latent',prefix),'w');
  makeTable(out,{'$\omega$','$\psi$','$E[\lambda]$'},q,fy);
  fclose(out);

  if (csv)
    writecsv(csv,{'Quantile','omega','psi','E[lambda]'}, ...
             [q(:,1) fy] );
  end
end

function makeTable(file,names,q,fy) 
  ind = round([0.1 0.25 0.5 0.75 0.9]*size(q,1))
  fprintf(file,' & %.2f ',[0.1 0.25 0.5 0.75 0.9]);
  fprintf(file,'\\\\ \\hline \n');
  for i=1:size(fy,2)
    fprintf(file,'%s',names{i})
    fprintf(file,'& %.3f ',fy(ind,i));
    fprintf(file,'\\\\ \n');
  end
end