function [hH, hB, hmu, history] = gfm_evaluate2(X, group, type, q, dropout, eps2, maxIter, omega, output, fast_version)
% The function to conduct the iterative algortihm
% Created by Wei Liu on 19/05/25.
% Copyright ? 2019 Wei Liu. All rights reserved.
% 
if(~exist('output', 'var'))
    output = 0;
end
if(~exist('fast_version', 'var'))
    fast_version = 0;
end
ind_set = unique(group);
ng = length(ind_set);
if ~isempty(setdiff(1:ng, ind_set))
    error('Id number of types must match type!')
end
if ng ~= size(type,1)
    error('The number of groups must match with length of type!');
end
hH = factorm(X, q);
gcell = cell(1, ng);
parfor j = 1:ng
    gcell{j} = find(group==j);
end
[n,p] = size(X);
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
    for j = 1:ng
        
        [B1] = localupdateB2(X, gcell{j}, hH, type(j,:));
         hhB = [hhB, B1];
    end
    hmu = hhB(1,:)';
    hB = hhB(2:end,:)';
    % ensure indentifiability.
    [B0tmp, ~] = qr(hB, 0);
    B0= B0tmp * diag(sort(sqrt(eig(hB'*hB)), 'descend'));
    
    sB = sign(B0(1,:));
    hB = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
    %hB(1:4,:), B(1:4,:)
    hBm = [hmu, hB];
    dB = norm(hBm - tmpB, 'fro')/norm(hBm, 'fro');
    tmpB = hBm;
    if output
        fprintf('-------------------------------------------\n')
        fprintf('---------- B updation is finished!---------\n')
    end
    % given B^(1), update H^(1)
    H4 = localupdateH2(X, gcell, hBm, type, dropout);
    if ng == 1 || fast_version == 1
        H5 = H4(:, 2:end);
    else
        H4(:,1) =1;
        H5 = localonestepH2(X, hBm, H4, gcell, type);
    end
    hH0 = H5;
    % H1(:,1:5)', H2(:,1:5)', H3(1:5,:), H(1:5,:), hH(1:5,:)
    [H0, ~] = qr(hH0, 0);
    hH1 = H0 * sqrt(n);
    sH0 = sign(hH0(1,:)).* sign(hH1(1,:));
    hH = hH1.* repmat(sH0,n,1);
    dH = norm(hH-tmpH, 'fro')/norm(hH, 'fro');
    tmpH = hH;
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