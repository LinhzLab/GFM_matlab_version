function [X, H, B, mu] = gendata(seed, n, p, type, q)
% type = {'homonorm', 'heternorm', 'pois', 'norm_pois', 'pois_bino', 'npb'}
if(~exist('n', 'var'))
    n = 300;
end
if(~exist('p', 'var'))
    p = 50;
end
if(~exist('type', 'var'))
    type = 'homonorm';
end
if(~exist('q', 'var'))
    q = 6;
end

switch type
    case 'homonorm'
        
        rng(1); % to fix B and mu
        ar_mat = 1*toeplitz(0.5.^(0:p-1)); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        % Zdecomp' * Zdecomp
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        % V = B'*B/p;
        mu = 0.4*randn(1,p);
        
        rng(seed);  % For reproducibility
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure F is unqiue (A2) condition
        
        X = repmat(mu, n, 1) + H * B' + mvnrnd(zeros(1,p), diag(ones(p,1)), n);
    case 'heternorm'
        rng(1);
        ar_mat = 1*toeplitz(0.5.^(0:p-1)); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        % Zdecomp' * Zdecomp
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1);
        mu = 2*randn(1,p);
        
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2); % ensure F is unqiue (A2) condition
        
        sigmas = 0.1+2*rand(p,1); % heteroskedasticity
        X = repmat(mu, n, 1) + H * B' + mvnrnd(zeros(1,p), diag(sigmas), n);
    case 'pois'
        rng(1); % since B is a fixed matrix.
        ar_mat = cov_mat(ones(p,1)*sqrt(4), 0.5); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);  % ensure B is unique.
        mu = 0.4*randn(1,p);
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        % H = mvnrnd(zeros(1,q),eye(q),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        % genarate X
        g2 = 1:p; % Poisson: exp
       
        mu2 = exp(H* B(g2,:)'+ repmat(mu(g2), n, 1)); % poisson distribution
        X = poissrnd(mu2);
    case 'norm_pois'
        rng(1); % since B is a fixed matrix.
        ar_mat = cov_mat(ones(p,1)*sqrt(4), 0.5); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);  % ensure B is unique.
        mu = 0.4*randn(1,p);
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        % H = mvnrnd(zeros(1,q),eye(q),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        % genarate X
        g1 = 1:p/2; % identity, normal
        g2 = (p/2+1):p; % Poisson: exp
        
        mu1 = H* B(g1,:)' + repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 1);
        mu2 = exp(H* B(g2,:)'+ repmat(mu(g2), n, 1)); % poisson distribution
        X2 = poissrnd(mu2);
        
        X = [X1 X2];
    case 'pois_bino'
        rng(1); % since B is a fixed matrix.
        ar_mat = cov_mat(ones(p,1)*sqrt(6), 0.5); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);  % ensure B is unique.
        mu = 0.4*randn(1,p);
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        
        % genarate X
        g1 = 1:floor(p/2); % poisson
        g2 = (floor(p/2)+1):p; % Poisson: exp
        
        
        mu1 = exp(H* B(g1,:)'+repmat(mu(g1), n, 1)); % poisson distribution
        X1 = poissrnd(mu1);
        
        mu2 = 1./(1 + exp(-H*B(g2,:)'- repmat(mu(g2), n, 1))); % binary distribution.
        
        X2 = binornd(1, mu2);
        X = [X1 X2];
    case 'npb'
        rng(1);
        Z = 1/3*randn(p,q);
        [B0tmp, ~] = qr(Z, 0);
        B0= B0tmp * diag(sort(sqrt(eig(Z'*Z)), 'descend'));
        
        sB = sign(B0(1,:));
        B = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
        mu = 0.4*randn(1,p); % Intercept
        % sB = sign(B(1,:));
        % B = B.* repmat(sB,p,1);  % ensure B is unique.
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        
        % genarate X
        g1 = 1:floor(p/3); % identity, normal
        g2 = (floor(p/3)+1):floor(2*p/3); % Poisson: exp
        g3 = (floor(2*p/3) + 1):p; % Bernoulli
        
        mu1 = H* B(g1,:)' + repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 1);
        
        mu2 = exp(H* B(g2,:)' + repmat(mu(g2), n, 1)); % poisson distribution
        X2 = poissrnd(mu2);
        
        mu3 = 1./(1 + exp(-H*B(g3,:)' - repmat(mu(g3), n, 1))); % binary distribution.
        X3 = binornd(1, mu3);
        
        X = [X1 X2 X3];
    case  'bino'
        q = 4;
        % generate H, B
        rng(1); % since B is a fixed matrix.
        ar_mat = cov_mat(ones(p,1)*sqrt(6), 0.5); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);  % ensure B is unique.
        mu = 0.4*randn(1,p); % Intercept
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        % genarate X
        g = 1:p;
        mu1 = 1./(1 + exp(-H*B(g,:)'- repmat(mu, n, 1))); % binary distribution.
        N = 1;
        X = binornd(N, mu1);
    case 'norm_bino'
        q = 8;
        % generate H, B
        rng(1); % since B is a fixed matrix.
        ar_mat = cov_mat(ones(p,1)*sqrt(4), 0.5); % p*AR(1) covariance matrix
        Z =  mvnrnd(zeros(1,p), ar_mat, n);
        [Zdecomp,~] = eig(Z*Z');
        B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1);  % ensure B is unique.
        mu = 0.4*randn(1,p); % Intercept
        % generate H, B
        rng(seed);
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        % H = mvnrnd(zeros(1,q),eye(q),n);
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue
        % genarate X
        g1 = 1:p/2; % identity, normal
        g2 = (p/2+1):p; % Poisson: exp
        
        mu1 = H* B(g1,:)'+ repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 1);
        mu2 = 1./(1 + exp(-H*B(g2,:)' - repmat(mu(g2), n, 1))); % binary distribution.
        X2 = binornd(1, mu2);
        X = [X1 X2];
end
end