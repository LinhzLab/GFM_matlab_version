function H2 = localonestepH2(X, Bm, Hm, gcell, type)
% function to conduct one-step updating for efficiently estimate latent factor matrix H.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
[n, p] = size(X);
B = Bm(:,2:end);
q = size(Hm, 2) - 1;
ng = size(type,1);

% mean matrix
mucell = cell(1,ng);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mucell{j} = Hm * Bm(gcell{j},:)';
    case 'poisson'
        mucell{j} = exp(Hm * Bm(gcell{j},:)');
    case 'binomial'
        Xj = X(:, gcell{j});
        ntrail_j = length(unique(Xj(:,1)))-1;
        mucell{j} =ntrail_j*1./(1 + exp(-Hm * Bm(gcell{j},:)'));
    end
end
% % Hessian matrix or information matrix
d2f = cell(n,1);
for i = 1:n
    Bng = zeros(q,q);
   for j = 1:ng
        switch type{j,1}
    case 'normal'
        %W = diag(1./ (std(X(:,gcell{j})).^2));
        Bng = Bng + B(gcell{j},:)'*B(gcell{j},:);
    case 'poisson'
        Bng = Bng + B(gcell{j},:)'* diag(mucell{j}(i,:)) * B(gcell{j},:);
        %Bng = Bng + (repmat(mucell{j}(i,:), q, 1)'.* B(gcell{j},:))'* B(gcell{j},:);
    case 'binomial'
        %Bng = Bng + (repmat(mucell{j}(i,:), q,1)' .* B(gcell{j},:))' *(repmat(1-mucell{j}(i,:), q,1)' .* B(gcell{j},:));
        Bng = Bng + (B(gcell{j},:))' * diag(mucell{j}(i,:).*(1-mucell{j}(i,:))) * B(gcell{j},:);

        end
   end
   d2f{i} = Bng;
end

% socre matrix
df2 = zeros(n, q);
for j = 1:ng
    df2 = df2 + (X(:, gcell{j})- mucell{j})* B(gcell{j},:);
end
cell1 = cell(n,1);

for i = 1:n
    cell1{i} = df2(i,:)';
end

H2update = cellfun(@(x, y) (x+ 1e-6*eye(q))\y, d2f, cell1, 'UniformOutput', 0);
% H2 = Hm - cell2array(H2update)';
H2updatemat = reshape(cell2mat(H2update), q,n);
H2 = Hm(:, 2:end) + H2updatemat';
