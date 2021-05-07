function obj = objfunc(H, B, X, omega, gcell, type)
% function to evaluate the value of the objective function (the conditional
% loglikelihood function).
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.

[n, p] = size(X);
eps1 = 1e-20;
Bh = H * B';
ng = size(type,1);
Q = zeros(size(X));
for j = 1:ng
    switch type{j,1}
    case 'normal'        
        Q(:, gcell{j}) = (X(:, gcell{j})- Bh(:, gcell{j})).^2;
    case 'poisson'
        me = exp(Bh(:,gcell{j}));
        Q(:, gcell{j}) = -(log(poisspdf(X(:, gcell{j}), me)+eps1));
    case 'binomial'
        me3 = 1 ./(1+exp(-Bh(:,gcell{j})));
        Q(:,gcell{j}) = -(X(:,gcell{j}).*log(me3+eps1) + (1-X(:,gcell{j})).* log(1-me3+eps1));
    end
end

obj = 1/n*omega*sum(sum(Q));
