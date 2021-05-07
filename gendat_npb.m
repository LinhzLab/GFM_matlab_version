function [X, H, B] = gendat_npb(i, n, p)
% This function generate simulated data with sample size n and dimension p.
% three exponential family distribution: Normal + Poisson + binary.
% Created by Wei Liu on 19/05/25.
% Copyright ? 2019 Wei Liu. All rights reserved.

if(~exist('p', 'var'))
   p = 50; 
end
if(~exist('n', 'var'))
   n = 300; 
end
q = 6;

% generate H, B
% Two methods to generate matrix B.
% Method 1.
% rng(1); % since B is a fixed matrix.
% ar_mat = cov_mat(ones(p,1)*sqrt(4), 0.5); % p*AR(1) covariance matrix
% Z =  mvnrnd(zeros(1,p), ar_mat, n);
% [Zdecomp,~] = eig(Z*Z');
% B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
% Method 2:
rng(1);
Z = 1/3*randn(p,q);
[B0tmp, ~] = qr(Z, 0);
B0= B0tmp * diag(sort(sqrt(eig(Z'*Z)), 'descend'));
    
sB = sign(B0(1,:));
B = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
% sB = sign(B(1,:));
% B = B.* repmat(sB,p,1);  % ensure B is unique.
% generate H, B
rng(i);
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g1 = 1:floor(p/3); % identity, normal
g2 = (floor(p/3)+1):floor(2*p/3); % Poisson: exp
g3 = (floor(2*p/3) + 1):p; % Bernoulli

mu1 = H* B(g1,:)'; % normal dstribution.
X1 = normrnd(mu1, 1); 
mu2 = exp(H* B(g2,:)'); % poisson distribution
X2 = poissrnd(mu2);
mu3 = 1./(1 + exp(-H*B(g3,:)')); % binary distribution.
X3 = binornd(1, mu3);

X = [X1 X2 X3];
end