%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Updated Date: Arpil 02, 2020
% Copyright (c) 2019, Liu Wei
% All rights reserved.
%%-------- CASE A. Including Intercepts
%% ------------ Example 1. Generalized factor model: all are normal with homoskedasticity
clear;
i = 4; p = 100; n = 100; q = 6; type = 'homonorm';
[X, H, B, mu] = gendata(i, n, p, type, q);

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; 
[hH, hB, ~, history] = gfm(X, group, type, q); % start to estimate.
%
q = []; % use PC criteria to dertermine q
[hH, hB, hmu, history] = gfm(X, group, type, q); % start to estimate.
[hH1, hB1] = factorm(X, 6, 0);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]


%% ------------ Example 2. Generalized factor model: all are normal with  heteroskedasticity
clear;
i = 4; p = 100; n = 100; q = 6;type = 'heternorm';
[X, H, B, mu] = gendata(i, n, p, type, q);

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; 
[hH, hB, ~, history] = gfm(X, group, type, q); % start to estimate.
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q); % start to estimate.

[hH1, hB1] = factorm(X, 6, 0);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 3. Generalized factor model: poisson + binary
clear;
i = 1; p = 100; n = 100; q=6; type = 'pois_bino';
[X, H, B, mu] = gendata(i, n, p, type, q);

[hH1, hB1] = factorm(X, 6);
hH1(1:6, 1:6)
H(1:6, 1:6)

% unified function
type = cell(2,2);
type{1,1} = 'poisson'; type{1,2} = 'log';
type{2,1} = 'binomial';  type{2,2} = 'probit'; % Also can choose 'logit'
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))]; % The first p/2 variables are of 'poisson'
% type, and the last p/2 variables are of 'binomial' type.
dropout=[2]; % Because 'binomial' type is weak signal relative to 'poisson' type, we add the 'dropout'
% step!
maxIter = 10;
q= 6; 
[hH, hB, hmu, history] = gfm(X, group, type, q, dropout, [], maxIter);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout); %%
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 4. Generalized factor model: normal + poisson
clear;
i =3; p = 100; n = 100;q=6; type = 'norm_pois';
[X, H, B, mu] = gendata(i, n, p, type, q);
B'*B
[hH1, hB1] = factorm(X, 6,0);
measurefun(H, hH1)
measurefun(B, hB1)
% unified function, well done, it is very very good.
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
type = {'normal', 'identity'; 'poisson', 'log'}; q= 6;
[hH, hB, ~, history] = gfm(X, group, type, q);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q); %%

%% ---------- Example 5. Generalized factor model: normal + poisson + binomial types
clear;
n = 200; p=200; i = 1;q=6; type = 'npb';
[X, H, B, mu] = gendata(i, n, p, type, q);
B'*B
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified functions test
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
q= 6; dropout=[3];omega = p^(-1);
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, [], [], omega);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout); %%
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% ---------- Example 6. Generalized factor model: pure binomial types

clear;
% generate data
p = 300; n = 300; i = 1;q=6; type = 'bino';
[X, H, B, mu] = gendata(i, n, p, type, q);
%----gfm unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'binomial';  
type{1,2} = 'logit'; 
q = 6; % q is given
[hH, hB, ~, history] = gfm(X, group, type, 6, 0, 1e-3, 6);
[hH1, hB1] = factorm(X, 6, 0);
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
% q is estimated; use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, [], 0, 1e-3, 6);
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% Fast algorithm for ultra-high-dimensional data
clear;
n = 5000; p=100000; i = 1;q=6; type = 'npb';
[X, H, B, mu] = gendata(i, n, p, type, q);
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified functions test
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
q= 6; dropout=[3];omega = p^(-1); output = 1; fast_version = 1;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, [], [], omega, [], output, fast_version);
history
% Compare the estimation performance by GFM and LFM.


%% -------- CASE B. Without Intercepts
%% ------------ Example 1. Generalized factor model: all are normal with homoskedasticity
clear;
i = 4; p = 100; n = 100; q = 6;
[X, B, H] = gendat_homonorm(i,n,p);

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; 
dropout = 0; dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% ------------ Example 2. Generalized factor model: all are normal with  heteroskedasticity
clear;
i = 4; p = 100; n = 100; q = 6;
[X, B, H] = gendat_heternorm(i,n,p);

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; 
dropout = 0; dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 3. Generalized factor model: poisson + binary
clear;
i = 1; p = 100; n = 100;
[X, H, B,hB] = gendat_pois_bino(i, n, p);
measurefun(B, hB)
[hH1, hB1] = factorm(X, 6);
hH1(1:6, 1:6)
H(1:6, 1:6)

% unified function
type = cell(2,2);
type{1,1} = 'poisson'; type{1,2} = 'log';
type{2,1} = 'binomial';  type{2,2} = 'probit'; % Also can choose 'logit'
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))]; % The first p/2 variables are of 'poisson'
% type, and the last p/2 variables are of 'binomial' type.
q = 6; dropout=[2]; % Because 'binomial' type is weak signal relative to 'poisson' type, we add the 'dropout'
dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 4. Generalized factor model: normal + poisson
clear;
i =3; p = 50; n = 50;
[X, H, B] = gendat_norm_pois(i, n, p);
B'*B
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified function, well done, it is very very good.
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
type = {'normal', 'identity'; 'poisson', 'log'}; q= 6;
dropout = 0; dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% ---------- Example 5. Generalized factor model: normal + poisson + binomial types
clear;
n = 200; p=200; i = 1;
[X, H, B] = gendat_npb(i, n, p);
B'*B
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified functions test
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
q= 6; dropout=[3];
 dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% ---------- Example 6. Generalized factor model: pure binomial types

clear;
% generate data
p = 300; n = 300; seed = 1;
[X, B, H] = gendat_binomial(seed,n,p);
%----gfm unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'binomial';  
type{1,2} = 'logit'; 
q = 6; % q is given
dropout = 0; dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%
q = []; % use PC criteria to dertermine q
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% Fast algorithm for ultra-high-dimensional data
clear;
n = 5000; p=100000; i = 1;
[X, H, B] = gendat_npb(i, n, p);
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified functions test
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
q= 6; dropout=[3];omega = p^(-1); output = 1; fast_version = 1;
dropout = 0; dc_eps = 1e-5;maxIter = 10; q_set = []; output = 0; fast_version=1; intercept = 0;
[hH, hB, ~, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)