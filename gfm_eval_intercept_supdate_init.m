function [hH, hB, hmu, history] = gfm_eval_intercept_supdate_init(X, group, type, q, dropout, eps2, maxIter, output)
% The function to conduct the iterative algortihm by simultaneously exert
% the identifiability conditions.
% Created by Wei Liu on 2020/01/18.
% Copyright ? 2019 Wei Liu. All rights reserved.
% 
if(~exist('output', 'var'))
    output = 0;
end
ind_set = unique(group);
ng = length(ind_set);
if ~isempty(setdiff(1:ng, ind_set))
    error('Id number of types must match type!')
end
if ng ~= size(type,1)
    error('The number of groups must match with length of type!');
end
hH = factorm(X, q, 0); % 0: specify X is uncentered.
gcell = cell(1, ng);

parfor j = 1:ng
    gcell{j} = find(group==j);
end
[n,p] = size(X);
omega = p^(-1);
% initialize
warning('off');
eps1 = 1e-4;
hB = 0;
dBm = Inf; dH = Inf; dOmega = Inf;dc =Inf;
dOmega = max([dBm, dH]);
tmpB = zeros(p,q+1);tmpH = hH; tmpc = 1e7;
% maxIter = 50;
k = 1;
tic;
while k <= maxIter && dOmega > eps1 && dc > eps2
    hhB = [];
    parfor j = 1:ng
        
        [B1] = localupdateB2(X, gcell{j}, hH, type(j,:));
        hhB = [hhB, B1];
    end
    hmu = hhB(1,:)'; % intercept estimate
    hB = hhB(2:end,:)'; % loading estimate.
    
    hBm = [hmu, hB];
    if output
        fprintf('-------------------------------------------\n')
        fprintf('---------- B updation is finished!---------\n')
    end
    % given B^(1), update H^(1)
    H4 = localupdateH2(X, gcell, hBm, type, dropout);
    hH0 = H4(:,2:end);
    hHBt = hH0 * hB';
    [U, S, V] = svd(hHBt, 'econ');
    hHs = sqrt(n) * U(:, 1:q); % sign not be identified
    hBs = V(:, 1:q) * S(1:q,1:q) / sqrt(n);
    sdhBs = diag(sign(hBs(1,:)) );
    hH = hHs * sdhBs;
    hB = hBs * sdhBs;
    dH = norm(hH-tmpH, 'fro') / norm(hH, 'fro');
    tmpH = hH;
    hBm = [hmu, hB];
    dB = norm(hBm - tmpB, 'fro')/norm(hBm, 'fro');
    tmpB = hBm;
    if output
        fprintf('-------------------------------------------\n')
        fprintf('---------- H updation is finished!---------\n')
        fprintf('-------------------------------------------\n')
    end 
    hHm = [ones(n,1), hH];
    dOmega = max([dB, dH]);
    c = objfunc(hHm, hBm, X, omega, gcell, type);
    dc = abs(c - tmpc)/abs(tmpc);
    tmpc = c;
    
    if output
        fprintf('Iter %d \n', k);
        fprintf('dB= %4f, dH= %4f,dc= %4f, c=%4f \n', dB, dH,dc, c);
    end
    history.dB(k) = dB; history.dH(k) = dH; history.dc(k)=dc; history.c(k)=c;history.realIter = k;
    k = k+1;
end
history.maxIter = maxIter;