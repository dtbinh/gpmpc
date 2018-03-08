

function [J, dJdp, dJdpStuct] = nvalue(p, m0, S0, dynmodel, policy, plant, cost, H, varargin)
%% Code

policy.p = rewrap(policy.p, p);            % overwrite policy.p with new parameters from minimize
p = unwrap(policy.p); dp = 0*p;
m = m0; S = S0; L = zeros(1,H);

if nargin > 8
    targetvec = varargin{1};
else
    tempvec = cost.z;
    targetvec = repmat(tempvec', H, 1);
end

if nargout <= 1                                       % no derivatives required
  
  for t = 1:H                                  % for all time steps in horizon
    cost.z = targetvec(t, :)';
    [m, S] = plant.prop(m, S, plant, dynmodel, policy);      % get next state
    L(t) = cost.gamma^t.*cost.fcn(cost, m, S);     % expected discounted cost
  end
else                                               % otherwise, get derivatives
  
  dmOdp = zeros([size(m0,1), length(p)]);
  dSOdp = zeros([size(m0,1)*size(m0,1), length(p)]);
  
  for t = 1:H                                  % for all time steps in horizon
    cost.z = targetvec(t, :)';
    [m, S, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
      plant.prop(m, S, plant, dynmodel, policy); % get next state
    
    dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;
    dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;
    
    [L(t), dLdm, dLdS] = cost.fcn(cost, m, S);              % predictive cost
    L(t) = cost.gamma^t*L(t);                                       % discount
    dp = dp + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )';
    
    dmOdp = dmdp; dSOdp = dSdp;                                 % bookkeeping
  end 
end

J = sum(L); 
dJdp= dp;
dJdpStuct = rewrap(policy.p, dp);
