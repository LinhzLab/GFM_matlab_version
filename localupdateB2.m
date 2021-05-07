function [B1] = localupdateB2(X, g1, hH, type1)
% function to update loading matrix B.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.

n = size(X, 1);
q = size(hH,2);
p1 = length(g1);
B1 = zeros(q+1, p1);
jg = 1:p1;
if strcmp(type1{1}, 'binomial')    
    parfor j = jg
        ntrail_j = length(unique(X(:,g1(j))))-1;
        B1(:,j) = glmfit(hH, [X(:,g1(j)), ntrail_j*ones(n,1)],type1{1},'link',type1{2}, 'constant', 'on');
    end
else
    parfor j = jg
    B1(:,j) = glmfit(hH, X(:,g1(j)),type1{1},'link',type1{2}, 'constant', 'on');
    end
end