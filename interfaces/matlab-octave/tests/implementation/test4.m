%%%
%
% verify implementation
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')
addpath('./lexlse_dual')

clear;clc

%% -----------------------------------------------------------------------

nVar = 10;
m    = [3,3,3];
r    = m;

nObj = length(m);

lexqr_struct = define_problem(nVar,m,r);

obj = lexqr_struct.obj;

%% append a terminal objective
if 0
    obj(nObj+1).A  = eye(nVar);
    obj(nObj+1).b  = zeros(nVar,1);
    obj(nObj+1).lb = obj(nObj+1).b;
    obj(nObj+1).ub = obj(nObj+1).b;
    nObj = nObj + 1;
end

%% -----------------------------------------------------------------------

options.regularization_type    = 7;
options.regularization_factors = ones(nObj,1);

[x,info,v0,as,d] = lexlsi(obj, options);

[xd,~,Xd] = lexlse_dual(obj, options.regularization_factors);

Xd - d.X_mu

for i = 1:nObj
    v{i}  = obj(i).A*x  - obj(i).b;
    vd{i} = obj(i).A*xd - obj(i).b;
    e(i)  = norm(v{i} - vd{i});
    e0(i) = norm(v0{i} - v{i});
end

fprintf('norm(v{i}-vd{i}) = ')
fprintf(' %e ', e)
fprintf('\n')

fprintf('norm(v{i}-v0{i}) = ')
fprintf(' %e ', e0)
fprintf('\n')

fprintf('norm(x-xd)       = %e\n',norm(x-xd))

%%%EOF
