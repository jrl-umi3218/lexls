%
% Copyright 2013-2021 INRIA
%

%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')

clear;clc

%% ===================================================================

nVar = 20;
m = [4,2,3,5];

nObj = length(m);

%% ===================================================================

for i=1:nObj
    obj(i).A  = randn(m(i),nVar);
    obj(i).lb = randn(m(i),1);
    obj(i).ub = obj(i).lb; % + rand(m(i),1);
end

%% ===================================================================

mu = 0.5*ones(nObj,1);

options.regularization_type    = 7;
options.regularization_factors = mu;

%% ===================================================================

[x_mu,info,v,as,d] = lexlsi(obj, options);

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

for k=1:nObj
    r{k} = obj(k).A*d.X_mu(:,k) - obj(k).ub;

    r{k} - v{k}
end

%%%EOF
