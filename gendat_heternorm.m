function [X, B, H] = gendat_heternorm(seed, n, p)
% This function generate simulated data with sample size n and dimension p.
% all are normal with heterogeneous vairances on j
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
if(~exist('p', 'var'))
   p = 50; 
end
if(~exist('n', 'var'))
   n = 300; 
end
q = 6;
rng(seed);  % For reproducibility
ar_mat = 1*toeplitz(0.5.^(0:p-1)); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
% Zdecomp' * Zdecomp
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1);
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2); % ensure F is unqiue (A2) condition
rng(1);
sigmas = 0.1+2*rand(p,1); % heteroskedasticity
X = H * B' + mvnrnd(zeros(1,p), diag(sigmas), n);
end