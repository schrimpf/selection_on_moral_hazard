
%% Reads in data and saves it in a nice structure
% file = name of data file
% drop = true if want to drop observations that change covariates 
% maxyear = use observations from maxyear and earlier
% lx = location of x's for psi, mu, and omega
% lxll = location of x's for kappa
% lxhs = location of x's for heteroskedasticity

%% data variables are
%1 id              long   %9.0g       id         Family ID
%2 year            int    %9.0g                  year
%3 choice          float  %9.0g                  
%4 avail1          float  %9.0g                  
%5 avail2          float  %9.0g                  
%6 avail3          float  %9.0g                  
%7 avail4          float  %9.0g                  
%8 avail5          float  %9.0g                  
%9 price1          float  %9.0g                  
%10 price2          float  %9.0g                  
%1 price3          float  %9.0g                  
%2 price4          float  %9.0g                  
%3 price5          float  %9.0g                  
%4 deduct1         float  %9.0g                  
%5 -deduct2         float  %9.0g                  
%6 -deduct3         float  %9.0g                  
%7 --m deduct4         float  %9.0g                  
%8 deduct5         float  %9.0g                  
%9 rate1           float  %9.0g                  
%20 rate2           float  %9.0g                  
%1 -rate3           float  %9.0g                  
%2 -rate4           float  %9.0g                  
%3 -rate5           float  %9.0g                  
%4 -max1            float  %9.0g                  
%5 -max2            float  %9.0g                  
%6 -max3            float  %9.0g                  
%7 -max4            float  %9.0g                  
%8 -max5            float  %9.0g                  
%9 -d_year_2004     byte   %8.0g                  year==2004
%30 -d_year_2005     byte   %8.0g                  year==2005
%1 -d_year_2006     byte   %8.0g                  year==2006
%2 -dcovg_tierr_2   byte   %8.0g                  covg_tierrecode==2
%3 -dcovg_tierr_3   byte   %8.0g                  covg_tierrecode==3
%4 -dcovg_tierr_4   byte   %8.0g                  covg_tierrecode==4
%5 -zgroup_2        byte   %8.0g                  group==2
%6 -zgroup_3        byte   %8.0g                  group==3
%7 -zgroup_4        byte   %8.0g                  group==4
%8 -zgroup_5        byte   %8.0g                  group==5
%9 -total           float  %9.0g                  Total Spending - Year
%40 age
%1  age2 (age^2/100)
%2  sex (1=female)
%3  tenure
%4  wage (income, in $1000)
%5  wage2 (income^2, in $10,000^2)
%6  family_size
%7 rs_fam
%8 rs_emp
%9 log(rs_fam)
%50 default
%1 rsfq1
%2 rsfq2
%3 rsfq3
%4 rsfq4

function data=loadData(file,drop,maxyear,lx,lxll,lxhs)
  if (nargin<2)
    drop = true;
  end
  if (nargin<3)
    maxyear = 2004;
  end

  %% read data from file
  rect = importdata(file,',');
  % only keep 2003 and 2004
  rect.data = rect.data(rect.data(:,2)<=maxyear,:);

  % a constant is automatically added to x's
  %lx{1} = [32:37 40 42 43 44 52 53 54]; % location of x's for omega (coverage tiers)
  %lx{2} = [32:37 40 42 43 44 52 53 54]; % location of x's for psi
  %lx{3} = [32:37 40 42 43 44 52 53 54]; % location of x's for lambda (allowed
                                  % to vary with t)
  avgX = zeros(size(lx{3}));
  avgX(end) = 0;
  %lxll = lx{3};
  %lxhs = lxll; % x's for heteroskedasticity
  lc = 3; % location of choice
  lav = 4:8; % location of avail
  lp = 9:13; % location of premiums
  ld = 14:18; % location of deductibles
  lr = 19:23; % location of copay rates
  lm = 24:28; % location of oop max
  lsp = 39; % total spending
  ldef = 50;

  data.id = rect.data(:,1);
  data.year = rect.data(:,2);
  uid = unique(data.id);
  data.N = numel(uid);
  uyear = unique(data.year);
  data.T = numel(uyear);
  %data.Z = rect.data(:,lz);
  prect = zeros(size(rect.data,2),data.T,data.N)*NaN;
  for t=1:data.T
    for i=1:data.N
      j = find(data.id==uid(i) & data.year==uyear(t));
      if (~isempty(j))
        prect(:,t,i) = rect.data(j,:)';
      end
    end
  end
  data.uid = mode(squeeze(prect(1,:,:)),1);
  data.choice = squeeze(prect(lc,:,:));
  data.totalSpend = squeeze(prect(lsp,:,:));
  data.default = squeeze(prect(ldef,:,:));
  avail = prect(lav,:,:);
  data.avail = avail;
  prem = prect(lp,:,:);
  data.prem = prem;
  deduct = prect(ld,:,:);
  data.deduct = deduct;
  copay = prect(lr,:,:);
  data.copay = copay;
  maxoop = prect(lm,:,:);
  data.maxoop = maxoop;
  
  % x's 
  for j=1:(numel(lx) - 1)
    xt = prect(lx{j},:,:);
    np = sum(~isnan(xt),2);
    xt(isnan(xt)) = 0;
    meanx = sum(xt,2)./np;
    data.x{j} = [ones(data.N,1) squeeze(meanx)'];
    % or mean(prect.data(lx{j},:,:),2)
  end
  xt = prect(lxll,:,:);
  np = sum(~isnan(xt),2);
  xt(isnan(xt)) = 0;
  meanx = sum(xt,2)./np;
  data.xll = [ones(data.N,1) squeeze(meanx)'];
  badx = ~all(squeeze(meanx(1:6,:,:))'==1 | squeeze(meanx(1:6,:,:))'==0,2);

  xt = prect(lxhs,:,:);
  np = sum(~isnan(xt),2);
  xt(isnan(xt)) = 0;
  meanx = sum(xt,2)./np;
  data.xhs = squeeze(meanx)';

  j=numel(lx);
  for t=1:data.T 
    xt = squeeze(prect(lx{j},t,:));
    dum = zeros(data.N,data.T-1);
    if (t>1) 
      dum(:,t-1) = 1;
    end
    data.x{j}{t} = [ones(data.N,1) xt' dum];
  end  
  % fill in missing x's with x's from other years -- these won't be used
  % in the estimates anyway, but NaN's cause problems 
  for t=1:data.T
    data.x{j}{1}(isnan(data.x{j}{1})) = ...
        data.x{j}{t}(isnan(data.x{j}{1}));
  end
  for t=1:data.T
    data.x{j}{t}(isnan(data.x{j}{t})) = ...
        data.x{j}{1}(isnan(data.x{j}{t}));
  end
  
  for t=1:data.T
    for j=2:size(data.xll,2)
      if (avgX(j-1))
        data.x{3}{t}(:,j) = data.xll(:,j);
      end
    end
  end

  if (sum(badx)>0 && drop)
    fprintf(['WARNING: %d observations change x\nDropping these observations.\n'], ...
            sum(badx))
    data=dropBad(data,badx);
  end    
  badx = zeros(size(data.x{1},1),1);
  for j = numel(data.x)-1
    badx = badx | any(isnan(data.x{j}),2); 
  end
  j = numel(data.x);
  for t = 1:numel(data.x{j})
    badx = badx | any(isnan(data.x{j}{t}),2);
  end
  badx = badx | any(isnan(data.xll),2);
  if (sum(badx)>0 && drop)
    fprintf('WARNING: %d observations have NaN for some x''s. Dropping them.\n', ...
            sum(badx))
    data=dropBad(data,badx);
  end    
  
  % identify and eliminate observations with impossible to rationalize
  % choice patterns (note now that mu_lambda changes over time, these choices
  % are not necessarily impossible to rationalize, so don't drop them.)
  %
  % bad = zeros(data.N,1);
  % for n=1:data.N
  %   for t1=1:(data.T-1)
  %     if ~isnan(data.choice(t1,n))
  %       for t2=(t1+1):data.T
  %         if (~isnan(data.choice(t2,n)) && data.choice(t1,n)~= ...
  %             data.choice(t2,n))
  %           bad(n) = isequalnan(data.avail(:,t1,n), data.avail(:,t2,n)) && ...
  %                    isequalnan(data.deduct(:,t1,n), data.deduct(:,t2,n)) && ...
  %                    isequalnan(data.prem(:,t1,n), data.prem(:,t2,n)) && ...
  %                    isequalnan(data.maxoop(:,t1,n), data.maxoop(:,t2,n));                                 
  %         end % if ~=
  %       end % for t2
  %     end % if ~isnan
  %   end % for t1
  %   if mod(n,100)==0
  %     fprintf('%6d of %6d\n',n,data.N);
  %   end
  % end % for i
  % if (sum(bad)>0 && drop)
  %   fprintf(['WARNING: %d observations have choices that cannot be ' ...
  %            'rationalized by model\nDropping these observations.\n'], ...
  %           sum(bad))
  %   data=dropBad(data,bad);
  % else
  %   data.dropped = [];
  % end
  %

end % function loadData

%% isequal treating NaN==NaN as true
function logical=isequalnan(x,y)
  logical = all(x==y | (isnan(x) & isnan(y)));
end
