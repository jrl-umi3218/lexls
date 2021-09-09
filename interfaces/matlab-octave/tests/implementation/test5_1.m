%
% Copyright 2013-2021 INRIA
%

%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%


addpath('/Users/drdv/Work/Installed/cpp/qpOASES-3.1.0/interfaces/matlab')
addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')
addpath('./lexlse_dual')

clear;clc

%% ===================================================================

BIG_NUMBER = 100;

k = 1;
obj(k).A  = 1;
obj(k).lb = -BIG_NUMBER;
obj(k).ub = 1;

k = 2;
obj(k).A  = 1;
obj(k).lb = 3;
obj(k).ub = BIG_NUMBER;

%% ===================================================================

nObj = length(obj);

mu = 1;

options.max_number_of_factorizations = 10;

options.regularization_type          = 7;
options.regularization_factors       = mu*ones(nObj,1);

%% ===================================================================


[x1,info,v1_,as,d] = lexlsi(obj, options);
[x2,~,v2_] = lex_sequence(obj,[],options.regularization_factors);
[x3,v3_]   = qpsequence(obj, mu);

[x1,x2,x3]

v1 = compute_violation(obj, x1);
v2 = compute_violation(obj, x2);
v3 = compute_violation(obj, x3);

fprintf('(status = %d) ', info.status);
fprintf(' %e ', [norm(v1{1}), norm(v2{1}), norm(v3{1})])
fprintf('\n')

tol = 1e-01;
if norm(v1{1}) < norm(v2{1}) - tol
    fprintf(' %e \n', abs(norm(v1{1}) - norm(v2{1})));

    %%fprintf(' %f ', e);
    %%fprintf('\n');
    %%keyboard
end

if abs(v2{1}-v3{1})>1e-12
    [norm(v2{1}), norm(v3{1})]
    keyboard
end

%%%EOF
