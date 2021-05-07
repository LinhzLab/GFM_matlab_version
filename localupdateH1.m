function [hH] = localupdateH1(X, gcell, hB, type, dropout)
% function to update latent factor matrix H.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.

n = size(X, 1);
q = size(hB, 2);
ng = size(type,1);
if all(dropout ~=0) && (~isempty(setdiff(dropout , 1:ng)))
    error('dropout setting is wrong!')
end
% dlcell = cell(ng, 2);
% for j = 1:ng
%     [dlcell{j,1}, dlcell{j,2}] = typematchlink(type{j});
% end
Harray = zeros(n,q, ng);
w = ones(1,n);
for j = 1:ng
    H2 = zeros(q,n);
    if strcmp(type{j,1}, 'normal')
        w = 1./ (std(X(:, gcell{j})).^2);
        parfor i = 1:n
            H2(:,i) = glmfit(hB(gcell{j},:), X(i,gcell{j})',type{j,1},'link',type{j,2}, 'constant', 'off', 'weights', w);
        
        end
    elseif strcmp(type{j,1}, 'binomial')
         p1 = length(gcell{j});
        parfor i = 1:n
            ntrail_i = length(unique(X(i,gcell{j})))-1;
           
            H2(:,i) = glmfit(hB(gcell{j},:), [X(i,gcell{j})',ntrail_i*ones(p1,1)],type{j,1},'link',type{j,2}, 'constant', 'off');
            
        end
    else
        parfor i = 1:n
            H2(:,i) = glmfit(hB(gcell{j},:), X(i,gcell{j})',type{j,1},'link',type{j,2}, 'constant', 'off');
            
        end
    end
    
    Harray(:,:,j) = H2';
end
idres = setdiff(1:ng, dropout);
hH = mean(Harray(:,:,idres), 3);