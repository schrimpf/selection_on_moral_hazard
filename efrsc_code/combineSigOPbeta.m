%% Combine parts of variance matrix into the variance matrix
function v = combineSigOPbeta(sigOP,gam,rho,T) 
  v = eye(2+T);
  v(1:2,1:2) = sigOP;
  v(3:(2+T),1) = (gam(1)*sigOP(1,1) + gam(2)*sigOP(1,2))*ones(T,1);
  v(1,3:(2+T)) =  v(3:(2+T),1)';
  v(3:(2+T),2) = (gam(1)*sigOP(1,2) + gam(2)*sigOP(2,2))*ones(T,1);
  v(2,3:(2+T)) =  v(3:(2+T),2)';
  rhoT = eye(T);
  for t=1:T
    for s=(t+1):T
      rhoT(t,s) =rho^(abs(t-s));
      rhoT(s,t) = rhoT(t,s);
    end
  end
  v(3:(2+T),3:(2+T)) = gam(1)^2*sigOP(1,1) + gam(2)^2*sigOP(2,2) + ...
      2*gam(1)*gam(2)*sigOP(1,2) + ...
      gam(3)^2/(1-rho^2) * rhoT;
end
