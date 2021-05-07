function Vr = ICriteria(X, hB, hH, r, group, type, criteria)
% PC or IC creteria to choose number of factors.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
if(~exist('criteria', 'var'))
    criteria = 'IC';
end

[n, p] = size(X);
omega = 1/p;
ind_set = unique(group);
ng = length(ind_set);
gcell = cell(1, ng);
for j = 1:ng
    gcell{j} = find(group==j);
end
c = objfunc(hH, hB, X, omega, gcell, type);
if strcmp(criteria, 'IC')
   Vr = [log(c) , r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)];
elseif strcmp(criteria, 'PC')
   Vr = [c , r*(n+p)/(n*p)*log(n*p/(n+p))];    
end
%Vr = [log(c) , r*(n+p)/(n*p)*log(min(sqrt(n), sqrt(p))^2)];
%Vr = [log(c) , r*(n+p)/(n*p)*log(n*p/(n+p))];
%Vr = [c , r*0.6*(n+p)/(n*p)*log(n*p/(n+p))];