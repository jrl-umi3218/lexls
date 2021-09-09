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

nObj = 2;
nVar = 10;
m = 3*ones(nObj,1);
r = m;

options.max_number_of_factorizations = 1000;

options.regularization_type    = 7;
options.regularization_factors = ones(nObj,1);

%% ===================================================================

for k=1:1

    if 1
	for i=1:nObj
	    if 0
		obj(i).A  = randn(m(i),nVar);
		obj(i).lb = randn(m(i),1);
		obj(i).ub = obj(i).lb + [rand(m(i)-1,1);0];
	    else
		if i == 1
		    obj(i).A  = randn(m(i),nVar);
		    obj(i).lb = randn(m(i),1);
		    obj(i).ub = obj(i).lb + [rand(m(i)-1,1);0];
		else
		    obj(i).A  = randn(m(i),nVar);
		    obj(i).lb = randn(m(i),1);
		    obj(i).ub = obj(i).lb;
		end
	    end
	end
    else
	load strange.mat
    end

    [x1,info,v1_,as1,d1] = lexlsi(obj, options);
    [x2,~,v2_] = lex_sequence(obj,[],options.regularization_factors);
    [x3,v3_]   = qpsequence(obj, 1);

    v1 = compute_violation(obj, x1);
    v2 = compute_violation(obj, x2);
    v3 = compute_violation(obj, x3);

    fprintf('(status = %d) ', info.status);
    fprintf(' %e ', [norm(v1{1}), norm(v2{1}), norm(v3{1})])
    fprintf('\n')

    for i = 1:nObj
	e(i) = norm(v1{i}) - norm(v2{i});
    end

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

end

%fprintf(' %e ', [norm(v1_{1}), norm(v2_{1}), norm(v3_{1})])
%fprintf('\n')

%%%EOF
