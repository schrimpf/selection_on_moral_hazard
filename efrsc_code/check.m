function ll=check(u,loglambda,lamlo,data) 
  % checks spending and latent variables. 
  lam = exp(loglambda)-(lamlo*ones(1,data.T))';
  tol = 1e-3;
  nbad = 0;
  ll = loglambda;
  for i=1:data.N
    om = exp(u(i,1));
    for t=1:data.T
      if (~isnan(data.totalSpend(t,i)) && data.choice(t,i)>0)
        s = data.totalSpend(t,i);
        d = data.deduct(data.choice(t,i),t,i);
        m = data.maxoop(data.choice(t,i),t,i);
        c = 0.1;
        if (s==0 && lam(t,i) > 0)
          fprintf('s==0 & lam>0 lamlo=%g\n',lamlo(i));
          ll(t,i) = log(lamlo(i)/2);
          nbad = nbad+1;
        end
        if (s>0 && s<d && norm(lam(t,i)-s)>tol)
          ll(t,i) = log(s+lamlo(i));
          fprintf('%d,%d lam~=s %g %g, need loglambda=%g\n',i,t,s,lam(t,i), ...
                  ll(t,i));
          nbad = nbad+1;    
        end
        if (s>d && d+c*(s-d)<m && norm(lam(t,i)+(1-c)*om-s)>tol)
          ll(t,i) = log(s+lamlo(i)-om*(1-c));
          fprintf('%d,%d lam+(1-c)*om~=s %g %g %g, need loglambda=%g\n', ...
                  i,t,s,lam(t,i)+(1-c)*om,u(i,1), ...
                  ll(t,i));
          nbad = nbad+1;    
        end
        if (d+c*(s-d)>m && norm(s-(lam(t,i)+om))>tol)
          ll(t,i) = log(s+lamlo(i)-om);
          fprintf('lam+om~=s %g %g, need loglambda=%g\n',s,lam(t,i)+om, ...
                  ll(t,i));
          nbad = nbad+1;    
        end
      end
    end
  end
  if (nbad>0)
    fprintf('%d bad\n',nbad);
    %pause();
  end
end