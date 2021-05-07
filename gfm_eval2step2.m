function [hH_final, hB_final, hmu_final, history] = gfm_eval2step2(X, group, type, q, dropout, eps2, maxIter, output,ICmethod)
% The function to compare the two methods to exert the identifiability conditions
% Created by Wei Liu on 2020/01/30.
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
try
hH_init = hH; hB_init= hB; hmu_init= hmu;
[hH_final, hB_final, hmu_final] = gfm_eval_intercept_osfinal2(X, hH_init, hB_init,hmu_init, group, type, ICmethod);

catch
    hH_final = hH; hB_final= hB; hmu_final= hmu;
end


end
