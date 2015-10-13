function lexqr_struct = lexqr_lambda(lexqr_struct)
%%%
%
% calls lexqr_lambda_obj.m to form the matrix of Lagrange multipliers
%

    %% -------------------------------------
    %% input handling
    %% -------------------------------------
    nObj = lexqr_struct.nObj;
    nCtr = lexqr_struct.nCtr;
    nVar = lexqr_struct.nVar;
    dim  = [lexqr_struct.obj_info(:).dim];

    %% -------------------------------------
    %% determine the right-hand-side
    %% -------------------------------------

    if isfield(lexqr_struct, 'x_mu')

	RT = lexqr_struct.Rr(:,1:end-1)';
	P  = lexqr_struct.P;
	rhs_all = inv(RT'*RT)*RT'*P'*lexqr_struct.x_mu;

    else
	rhs_all = zeros(nVar,1);
    end

    %% -------------------------------------

    L = zeros(nCtr,nObj);
    for k = 1:nObj
	ind = 1:sum(dim(1:k));
	L(ind,k) = lexqr_lambda_obj(lexqr_struct,k,rhs_all);
    end

    %% DONT FORGET TO SCALE BY mu(k)

    lexqr_struct.lambda = L;

%%%EOF
