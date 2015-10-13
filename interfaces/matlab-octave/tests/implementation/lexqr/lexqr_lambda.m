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
    dim  = [lexqr_struct.obj_info(:).dim];
    %% -------------------------------------

    L = zeros(nCtr,nObj);
    for k = 1:nObj
	ind = 1:sum(dim(1:k));
	L(ind,k) = lexqr_lambda_obj(lexqr_struct,k);
    end

    lexqr_struct.lambda = L;

%%%EOF
