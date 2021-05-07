function [X, H, B,hB] = gendat_pois_bino(i, n, p)
% This function generate simulated data with sample size n and dimension p.
% two exponential family distribution: Poisson and binary.
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
% generate H, B
rng(1); % since B is a fixed matrix.
ar_mat = cov_mat(ones(p,1)*sqrt(6), 0.5); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);  % ensure B is unique.
% generate H, B
rng(i);
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g1 = 1:floor(p/2); % poisson
g2 = (floor(p/2)+1):p; % Poisson: exp

%mu1 = H* B(g1,:)'; % normal dstribution.
%X1 = normrnd(mu1, 0.01); 
mu1 = exp(H* B(g1,:)'); % poisson distribution
X1 = poissrnd(mu1);
mu2 = 1./(1 + exp(-H*B(g2,:)')); % binary distribution.
X2 = binornd(1, mu2);
X = [X1 X2];
B1 = [];
for j = g1
    b1 = glmfit(H, X(:,j),'poisson','link','log', 'constant', 'off');
%b1 = glmfit(H, X(:,j),'normal', 'constant', 'off');
B1 = [B1, b1];
end
j = 0;
B2 = [];
for j =g2
%b2 = glmfit(H, X(:,j),'poisson','link','log', 'constant', 'off');
b2 = glmfit(H, X(:,j), 'binomial', 'link', 'logit', 'constant', 'off');
B2 = [B2, b2];
end
hB = [B1,B2]';
end