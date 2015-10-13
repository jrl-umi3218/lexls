%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')

clear;clc

%% -----------------------------------------------------------------------

nObj = 4;
nVar = 15;
m = 5*ones(nObj,1);
r = m-2;

options.regularization_type    = 1;
options.regularization_factors = 0.5*ones(nObj,1);

%% ===================================================================

lexqr_struct = define_problem(nVar,m,r);

lexqr_struct.options = options;

%% ===================================================================

lexqr_struct = lexqr_form(lexqr_struct);
lexqr_struct = lexqr_solve(lexqr_struct);
lexqr_struct = lexqr_residual(lexqr_struct);
lexqr_struct = lexqr_explicit(lexqr_struct);
lexqr_struct = lexqr_compact(lexqr_struct);

[x_mu,info,v,as,d] = lexlsi(lexqr_struct.obj, options);
[~,~,~,~,d1] = lexlsi(lexqr_struct.obj);

%% When I use regularization in lexlse, the computation of Lambda changes. This is because the
%% residuals change. But I might be doing something wrong with the singular part of the RHS vector -
%% to check.

%lexqr_struct.x_mu = x_mu;

lexqr_struct = lexqr_lambda(lexqr_struct);

%% ===================================================================

L0 = [];
for k=1:lexqr_struct.nObj
    L0 = [L0;d.lambda{k}];
end

L1 = [];
for k=1:lexqr_struct.nObj
    L1 = [L1;d1.lambda{k}];
end

%lexqr_struct.lambda - L0

L1 - L0

%%%EOF
