function [B1] = localupdateB1(X, g1, hH, type)
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
B1 = zeros(q, p1);
jg = 1:p1;
if strcmp(type{1}, 'binomial')    
    parfor j = jg
        ntrail_j = length(unique(X(:,g1(j))))-1;
        B1(:,j) = glmfit(hH, [X(:,g1(j)), ntrail_j*ones(n,1)],type{1},'link',type{2}, 'constant', 'off');
    end
else
    parfor j = jg
    B1(:,j) = glmfit(hH, X(:,g1(j)),type{1},'link',type{2}, 'constant', 'off');
    end
end