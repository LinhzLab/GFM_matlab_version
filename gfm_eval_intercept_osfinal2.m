function [hH_final, hB_final, hmu_final, objvalue] = gfm_eval_intercept_osfinal2(X, hH_init, hB_init,hmu_init, group, type, ICmethod)
% The function to conduct the iterative algortihm
% osfinal: one-step to final estimates
% It is used to validate the effeciency gain by one-step update.
% Created by Wei Liu on 2020/01/30.
% Updated date: 2020-01-30
% Copyright ? 2020 Wei Liu. All rights reserved.
%
if(~exist('ICmethod', 'var'))
    ICmethod = 'SVD';
end
n = size(hH_init,1);
p = size(hB_init,1);
omega = p^(-1);
ind_set = unique(group);
ng = length(ind_set);
if ~isempty(setdiff(1:ng, ind_set))
    error('Id number of types must match type!')
end
if ng ~= size(type,1)
    error('The number of groups must match with length of type!');
end
if ng == 1 % If there is only one type, the initial estimator is the efficient estimator.
    hH_final = hH_init;
    hB_final = hB_init;
    hmu_final = hmu_init;
    gcell{1} = 1:p;
    objvalue = objfunc([ones(n,1), hH_final], [hmu_final, hB_final], X, omega, gcell, type);
else
    gcell = cell(1, ng);
    for j = 1:ng
        gcell{j} = find(group==j);
    end
    q = size(hH_init,2);
    % update H
    hBm_init = [hmu_init, hB_init];
    hHm_init = [ones(n,1), hH_init];
    H5 = localonestepH2(X, hBm_init, hHm_init, gcell, type);
  
    % Bm2 = localonestepB(X, hBm_init, hHm_init, gcell, type);
    % Bm_init = hBm_init; Hm_init= hHm_init;
%     hmu = Bm2(:,1);
%     hB = Bm2(:,2:end);
    hhB = [];
    for j = 1:ng
        
        % [B1] = localupdateB2(X, gcell{j}, hH_init, type(j,:));
        [B1] = localupdateB2(X, gcell{j}, H5, type(j,:));
        hhB = [hhB, B1];
    end
    hmu = hhB(1,:)';
    hB = hhB(2:end,:)';
    hBm_init = [hmu, hB];
    hHm_init = [ones(n,1), H5];
    H5 = localonestepH2(X, hBm_init, hHm_init, gcell, type);
    switch ICmethod
        case 'SVD'
            hHBt = H5 * hB';
            [U, S, V] = svd(hHBt, 'econ');
            hHs = sqrt(n) * U(:, 1:q); % sign not be identified
            hBs = V(:, 1:q) * S(1:q,1:q) / sqrt(n);
            sdhBs = diag(sign(hBs(1,:)) );
            hH_final = hHs * sdhBs;
            hB_final = hBs * sdhBs;
            hmu_final = hmu;
        case 'sep'
            hB = hhB(2:end,:)';
            [B0tmp, ~] = qr(hB, 0);
            B0= B0tmp * diag(sort(sqrt(eig(hB'*hB)), 'descend'));
            sBMat = diag(sign(B0(1,:)));
            hB_final = B0 * sBMat;
            [H0, ~] = qr(H5, 0);
            hH1 = H0 * sqrt(n);
            sH0 = sign(H5(1,:)).* sign(hH1(1,:));
            hH_final = hH1.* repmat(sH0,n,1);
            hmu_final = hmu;
    end
    objvalue = objfunc([ones(n,1), hH_final], [hmu_final, hB_final], X, omega, gcell, type);
end

end