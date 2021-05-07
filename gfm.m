function [hH, hB, hmu, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, q_set, output, fast_version, intercept)
%% This function is used to conduct the Generalized Factor Model.
% Author: Liu Wei(weiliu@smail.swufe.edu.cn)
% Version: gmf-0.2
% Created Date: 2018-10-22
% Updated Date: 2019-05-24
% Iterative algorithm to estimate the factors H and loading matrix B in
% generalized factor model. This is conducted
% parallelly across variables and individuals.
% 
% Copyright ? 2019 Wei Liu. All rights reserved.
%--------INPUT:
% --X: a matrix with dimension of n*p(p=(p1+p2+..+p_d)), observational mixed data
% matrix, d is the types of variables, p_j is the dimension of j-th type
% variable.
% --group: a vector with length equal to p, specify each column of X belonging to which group.
% --type: a d-by-2 cell variable, specify the type of variable and link function on each
% group. For example, type={'poisson', 'log'; 'binomial', 'probit'}, and it
% is referred to the help file of glmfit() function for more details.
% --q: a positive integer or empty, specify the number of factors. If q is
% empty, then IC criteria is used to dertemined q automatically.
% --dropout: a proper subset of [1, 2, ..., d],  specify which group to be
% dropped in obtaining the initial estimate of factor matrix H, and the aim
% is to ensure the convergence of algorithm leaded by weak signal variable
% types. Optional parameter with default as 0, no group dropping.
% -- maxIter: a positive integer, specify the times of iteration. Optional
% parameter with default as 50.
% -- dc_eps: a positive real number, specify the tolerance of varing quantity of objective
% function in the algorithm. Optional parameter with default as 1e-6.
% -- q_set: a positive integer set, specify the candidates of factor number
% q, (optional) default as [1:8] according Bai,2013.
% -- output: a logical value with true or false, specify whether ouput the
% mediate information in iteration process, (optional) default as false.
% -- fast_version: a integer value with 1 or 0, fast_version = 1: use the fast
% algorithm which omit the one-step updating, but it cannot ensure the estimated
% efficieny; fast_version = 0: use the original algorithm; (optional)
% default as 0; 
% -- intercept: a integer value with 1 or 0, indicates whether estimate the
% intercept term in GFM, default as 1.
%-------OUTPUT:
%--hH: a n*q matrix, the estimated factor matrix.
%--hB: a p*q matrix, the estimated loading matrix.
%--history: a structure variable including the following 9 fileds: 
% dB: the varied quantity of B in each iteration;
% dH: the varied quantity of H in each iteration;
% dc: the varied quantity of the objective function in each iteration; 
% c: the objective value in each iteration; 
% realIter: the real iterations to converge; 
% maxIter: the tolerance of maximum iterations;
% dc_eps: the tolerance of the varied quantity of objective fucntion;
% dropout: which group is dropped out, 0 indicates no group dropped out;
% q: the used or estimated factor number.
%-------EXAMPLES
% see gfm_example.m file for examples of gfm() function.
%----------------------------------------------------------------------------
if(~exist('dropout', 'var') || isempty(dropout))
    dropout=0;
end
if(~exist('dc_eps', 'var')|| isempty(dc_eps))
    dc_eps = 1e-6; 
end
if(~exist('maxIter', 'var')|| isempty(maxIter))
    maxIter = 50; 
end
X = double(X);
[n, p] = size(X);
omega = p^(-1/2);


if(~exist('fast_version', 'var'))
    fast_version = 0;
end
if(~exist('intercept', 'var')|| isempty(intercept))
    intercept = 1; 
end
if(~exist('output', 'var') || isempty(output))
    output = false;
end
if(~exist('q', 'var') || isempty(q))
  
    if(~exist('q_set', 'var') || isempty(q_set))
       q_set = 5:7; % Here, just take {5, 6, 7} as an example.
    end
    fprintf('Start to determine the factor number q ....\n')
    q = singleIC(X, group,type, q_set, dropout, dc_eps, 10, omega, output, fast_version, intercept);
 
    fprintf('The factor number q is estimated as %d . \n', q);
end


fprintf('Starting the iterative algorithm ....\n')
if intercept == 0
[hH, hB, history] = gfm_evaluate(X, group, type, q, dropout, dc_eps, maxIter, omega, output, fast_version);
 hmu = zeros(p,1);
elseif intercept == 1
    if fast_version == 1
        [hH, hB, hmu, history] = gfm_eval_intercept_init(X, group, type, q, dropout, dc_eps, maxIter,output);
    else
         [hH, hB, hmu, history] = gfm_eval2step2(X, group, type, q, dropout, dc_eps, maxIter, output);
    end
  % [hH, hB, hmu, history] = gfm_evaluate2(X, group, type, q, dropout, dc_eps, maxIter, omega, output, fast_version);
end
fprintf('Finishing the iterative algorithm ....\n')
history.dc_eps = dc_eps; history.dropout = dropout; history.q = q;
end
