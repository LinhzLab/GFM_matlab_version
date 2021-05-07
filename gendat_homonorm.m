function [X, B, H] = gendat_homonorm(seed,n,p)
%% linear factor model to generate data.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
q = 6;
if(~exist('n', 'var'))
    n = 300;
end
if(~exist('p', 'var'))
    p = 50;
end
%-------------------- Case I -----------------------%
rng(seed);  % For reproducibility
ar_mat = 1*toeplitz(0.5.^(0:p-1)); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
% Zdecomp' * Zdecomp
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
% V = B'*B/p;

H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure F is unqiue (A2) condition

X = H * B' + mvnrnd(zeros(1,p), diag(ones(p,1)), n);


