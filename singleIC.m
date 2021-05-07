function q = singleIC(X_resamp, group,type, q_set, dropout, eps2, maxIter, omega, output, fast_version,intercept)
% function to choose factor number from a candidate set by PC or IC criteria.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.


%dropout = 0; dc_eps = 1e-6;maxIter=50;
if(~exist('q_set', 'var') || isempty(q_set))
   q_set = 1:10;
end
n = size(X_resamp, 1);
q_num = length(q_set);
allhH = cell(q_num,1);
allhB = cell(q_num,1);
if intercept == 0
for r = 1:q_num
    %fprintf('r = %d \n ', r);
    [hH, hB] = gfm_evaluate(X_resamp, group, type, q_set(r), dropout, eps2, maxIter, omega, output, fast_version);
    allhH{r} = hH;
    allhB{r} = hB;
end
elseif intercept == 1
        for r = 1:q_num
    %fprintf('r = %d \n ', r);
    [hH, hB, hmu] = gfm_eval2step(X_resamp, group, type, q_set(r), dropout, eps2, maxIter, output);
    %[hH, hB, hmu] = gfm_evaluate2(X_resamp, group, type, q_set(r), dropout, eps2, maxIter, omega, output, fast_version);
    allhH{r} = [ones(n, 1), hH];
    allhB{r} = [hmu, hB];
        end
end
Vr = zeros(2, q_num);
for r = 1:q_num
    %fprintf('r = %d \n ', r)
    Vr(:, r) = ICriteria(X_resamp, allhB{r}, allhH{r}, r, group, type);
end
IC = sum(Vr);
[~, id] = min(IC);
q = q_set(id);