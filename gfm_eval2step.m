function [hH_final, hB_final, hmu_final, history] = gfm_eval2step(X, group, type, q, dropout, eps2, maxIter, output,ICmethod)
% The function to conduct the iterative algortihm
% Created by Wei Liu on 2020/01/28.
% Updating date: 2020-01-30
% Copyright ? 2019 Wei Liu. All rights reserved.
% ICmethod: method to achieve the identifiability condition: there are
% "sep" and "SVD"
if(~exist('output', 'var'))
    output = 0;
end
if(~exist('ICmethod', 'var'))
    ICmethod = 'sep'; % seperately to achieve the identifiability
end
switch ICmethod
    case 'sep'
       [hH, hB, hmu, history] = gfm_eval_intercept_init(X, group, type, q, dropout, eps2, maxIter,output);
    case 'SVD'
        [hH, hB, hmu, history] = gfm_eval_intercept_supdate_init(X, group, type, q, dropout, eps2, maxIter,output);

end
 
hH_init = hH; hB_init= hB; hmu_init= hmu;
[hH_final, hB_final, hmu_final] = gfm_eval_intercept_osfinal(X, hH_init, hB_init,hmu_init, group, type);
end
