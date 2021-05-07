function [X, B, H] = gendat_binomial(seed,n,p)
% This function generate simulated data with sample size n and dimension p.
% pure Binomial distribution
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
if(~exist('n', 'var'))
    n = 300;
end
if(~exist('p', 'var'))
    p = 50;
end
%-------------------- Case I -----------------------%

q = 6;
% generate H, B
rng(1); % since B is a fixed matrix.
ar_mat = cov_mat(ones(p,1)*sqrt(6), 0.5); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);  % ensure B is unique.
% generate H, B
rng(seed);  % For reproducibility
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g = 1:p;
mu = 1./(1 + exp(-H*B(g,:)')); % binary distribution.
N = 2;
X = binornd(N, mu);

