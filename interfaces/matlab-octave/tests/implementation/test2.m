%%%
%
%
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')
addpath('./lexlse_dual')

clear;clc

nObj = 7;
nVar = 20;
m    = 5*ones(nObj,1);
r    = m-2;
mu   = 0.5*ones(nObj,1);

nCtr = sum(m);

lexqr_struct = define_problem(nVar,m,r);
obj = lexqr_struct.obj;

[Ab, m] = obj2array(obj);
A = Ab(:,1:end-1);

%% ---------------------------------------------------------
options.regularization_type     = 1;
options.regularization_factors  = mu;
%% ---------------------------------------------------------

[x0,info,v,as,d] = lexlsi(obj, options);

L0 = [];
for i=1:nObj
    L0 = [L0; d.lambda{i}];
end

%% ---------------------------------------------------------

[xd, Ld, Xd] = lexlse_dual(obj,mu);

muXd = [];
for i=1:nObj
    muXd = [muXd, mu(i)^2*Xd(:,i)];
end

A'*Ld + muXd

[x0,xd,x0-xd]

%% ---------------------------------------------------------

%%%EOF
