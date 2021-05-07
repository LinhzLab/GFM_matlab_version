%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         Jan. 13, 2018
% Copyright (c) 2018, Liu Wei
% All rights reserved.
% ----------------------------------------------------------------------------------
% This function generate two type of correlation matrix.
% Syntax: covmat = cov_mat(sdvec, rho, type)
% INPUT ARGUMENTS:
% -sdvec:--a positive vector, standard deviation of each random variable.
% -rho: --a value between 0 and 1, a baseline vlaue of correlation coefficient.
% -type:--a character, specify the type of correlation matrix and only include
%        'toeplitz' and 'identity' in current version, default as
%        'toeplitz'(AR(1) correlation matrix).
% OUTPUT ARGUMENTS:
% -covmat:-- a covariance matrix with a type of specified structure.
function covmat = cov_mat(sdvec, rho, type)
    if(~exist('type', 'var'))
        type = 'toeplitz';
    end
    p = length(sdvec);
    covmat = zeros(p,p);
    if strcmp(type, 'toeplitz')
     cormat = toeplitz(rho.^(0:p-1));
    end
    if strcmp(type, 'identity')
     cormat = eye(p);
     cormat(cormat==0) = rho;
    end
    for i= 1:p
        for j = 1:p
            covmat(i,j) = sdvec(i)*sdvec(j)*cormat(i,j);
        end
    end
end