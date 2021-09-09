%
% Copyright 2013-2021 INRIA
%

%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')
addpath('./lexlse_dual')

clear;clc

%% -----------------------------------------------------------------------

if 0
    nObj = 4;
    nVar = 25;
    m = 5*ones(nObj,1);
    r = m-2;
else
    nVar = 9;
    m = [4,2,3,5];
    r = m-2;

    nObj = length(m);
end

options.regularization_type    = 7;
options.regularization_factors = 0.5*ones(nObj,1);

%% ===================================================================

lexqr_struct = define_problem(nVar,m,r);
obj = lexqr_struct.obj;
mu  = options.regularization_factors;
Ab  = obj2array(obj);
A   = Ab(:,1:end-1);

lexqr_struct.options = options;

%% ===================================================================

lexqr_struct = lexqr_form(lexqr_struct);
lexqr_struct = lexqr_solve(lexqr_struct);
lexqr_struct = lexqr_residual(lexqr_struct);
lexqr_struct = lexqr_explicit(lexqr_struct);
lexqr_struct = lexqr_compact(lexqr_struct);

[x_mu,info,v,as,d] = lexlsi(lexqr_struct.obj, options);

lexqr_struct.X_mu = d.X_mu;

csm = [0;cumsum(m(:))];
for i=1:nObj
    ind = csm(i)+1:csm(i+1);
    lexqr_struct.residual_mu{i} = d.residual_mu(ind);
end

[xd, Ld, Xd] = lexlse_dual(obj,mu);

lexqr_struct.Xd = Xd;

lexqr_struct = lexqr_lambda(lexqr_struct);

%% ===================================================================

L0 = [];
for i=1:nObj
    L0 = [L0;d.lambda{i}];
end

if 1
    muXd = [];
    muLambda = [];
    for i=1:nObj
	muXd = [muXd, mu(i)^2*Xd(:,i)];
    end

    fprintf('norm(A''*Ld + muXd)             = %e\n', norm(A'*Ld + muXd))
    fprintf('norm(A''*L0 + muXd)             = %e\n', norm(A'*L0 + muXd))
    fprintf('norm(Xd - d.X_mu)              = %e\n', norm(Xd - d.X_mu))
    fprintf('norm(lexqr_struct.lambda - L0) = %e\n', norm(lexqr_struct.lambda - L0))
end

%%%EOF
