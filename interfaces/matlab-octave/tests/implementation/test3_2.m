%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%

addpath('/Users/drdv/Work/Installed/cpp/qpOASES-3.1.0/interfaces/matlab')
addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath('./lexlse_dual')

clear;clc

%% ===================================================================

if 1
    k = 1;
    obj(k).A  = 1;
    obj(k).lb = 1;
    obj(k).ub = 1;
    obj(k).b  = obj(k).ub;

    k = 2;
    obj(k).A  = 1;
    obj(k).lb = 3;
    obj(k).ub = 3;
    obj(k).b  = obj(k).ub;
else
    k = 1;
    obj(k).A  = 1;
    %%obj(k).lb = 1;
    %%obj(k).ub = 2;
    obj(k).lb = -100;
    obj(k).ub = 1;
end

%% ===================================================================

nObj = length(obj);
mu = 1*ones(nObj,1);

%%options.max_number_of_factorizations = 200;

options.regularization_type    = 7;
options.regularization_factors = mu;

options.output_file_name = 'debug1.m';

%% ===================================================================

[x1,info,v1,as,d] = lexlsi(obj, options);

[x2,~,v2] = lex_sequence(obj,[],options.regularization_factors);

L0 = [];
A  = [];
for i=1:nObj
    L0 = [L0;d.lambda{i}];
    A = [A;obj(i).A];
end

muX_mu = [];
for i=1:nObj
    muX_mu = [muX_mu, mu(i)^2*d.X_mu(:,i)];
end

fprintf('norm(A''*L0 + muX_mu) = %e\n', norm(A'*L0 + muX_mu))

[xd, Ld, Xd] = lexlse_dual(obj,mu);

%% ===================================================================
%% CVX (testing with 3 objectives)
%% ===================================================================

if 0
    cvx_precision best
    cvx_clear

    cvx_begin quiet
    variables x1(1) r1(1)
    dual variables L11

    minimize(0.5*r1'*r1 + 0.5*mu(1)^2*x1'*x1)
    subject to
    L11: r1 == obj(1).A*x1 - obj(1).b
    cvx_end

    fprintf('primal(1): cvx_status = %s\n',cvx_status);

    cvx_begin quiet
    variables x2(1) r2(1)
    dual variables L21 L22

    minimize(0.5*r2'*r2 + 0.5*mu(2)^2*x2'*x2)
    subject to
    L21: r1 == obj(1).A*x2 - obj(1).b
    L22: r2 == obj(2).A*x2 - obj(2).b
    cvx_end

    fprintf('primal(2): cvx_status = %s\n',cvx_status);

    L1 = [L11 L21;0 L22];

    L0
    L1
end

%%%EOF
