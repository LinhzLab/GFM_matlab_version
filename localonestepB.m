function Bm2 = localonestepB(X, Bm_init, Hm_init, gcell, type)
% function to conduct one-step updating for efficiently estimate latent factor matrix H.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
[~, p] = size(X);
B = Bm_init(:,2:end);
q = size(Hm_init, 2) - 1;
ng = size(type,1);

% mean matrix
Bm2 = zeros(p, q+1);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mutypej = Hm_init * Bm_init(gcell{j},:)'; % n*p1
        scorej = Hm_init' * (X(:, gcell{j}) -mutypej); % (q+1) * p1
        Hesj = - Hm_init'* Hm_init;
        Bm2(gcell{j},:) = Bm_init(gcell{j},:) + scorej' / Hesj;
    case 'poisson'
        jvec = gcell{j};
        mutypej = exp(Hm_init * Bm_init(jvec,:)');
        scorej = Hm_init' * (X(:, jvec) -mutypej); % (q+1) * p1
        p1 = length(jvec);
        for jl = 1:p1
            Hestypejjl = - Hm_init'* diag(mutypej(:,jl)) * Hm_init;
            Bm2(jvec(jl),:) = Bm_init(jvec(jl),:) + scorej(:,j)' / Hestypejjl;
        end
        
    case 'binomial'
        jvec = gcell{j};
        Xj = X(:, jvec);
        ntrail_j = length(unique(Xj(:,1)))-1;
        mutypej =ntrail_j*1./(1 + exp(-Hm_init * Bm_init(jvec,:)'));
        scorej = Hm_init' * (Xj -mutypej); % (q+1) * p1
        p1 = length(jvec);
        
        for jl = 1:p1
            Hestypejjl = - Hm_init'* diag(mutypej(:,jl).*(1-mutypej(:,jl)))  * Hm_init;
            Bm2(jvec(jl),:) = Bm_init(jvec(jl),:) + scorej(:,j)' / Hestypejjl;
        end
    end
end

