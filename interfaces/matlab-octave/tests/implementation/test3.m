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

if 1
    nObj = 4;
    nVar = 25;
    m = 5*ones(nObj,1);
    r = m-2;
else
    nObj = 3;
    nVar = 15;
    m = [2,3,4];
    r = m;
end

options.regularization_type    = 7;
options.regularization_factors = 1*ones(nObj,1);

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

%% When I use regularization in lexlse, the computation of Lambda changes. This is because the
%% residuals change. But I might be doing something wrong with the singular part of the RHS vector -
%% to check. HOW TO COMPUTE THE RESIDUAL?

lexqr_struct.X_mu = d.X_mu;
csm = [0;cumsum(m(:))];
for i=1:nObj
    ind = csm(i)+1:csm(i+1);
    lexqr_struct.residual_mu{i} = d.residual_mu(ind);
end

[xd, Ld, Xd] = lexlse_dual(obj,mu); %% use ^2?

lexqr_struct = lexqr_lambda(lexqr_struct);

%% ===================================================================

if 1
    muXd = [];
    muLambda = [];
    for i=1:nObj
	muXd = [muXd, mu(i)^2*Xd(:,i)];
	muLambda = [muLambda, mu(i)^2*lexqr_struct.lambda(:,i)];
    end

    fprintf('norm(A''*Ld + muXd)       = %e\n', norm(A'*Ld + muXd))
    fprintf('norm(A''*muLambda + muXd) = %e\n', norm(A'*muLambda + muXd))
    fprintf('norm(Xd - d.X_mu)        = %e\n', norm(Xd - d.X_mu))

end

lexqr_struct.lambda - Ld %% could be different

%%%EOF
