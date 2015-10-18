function Lambda = lexqr_lambda_obj(lexqr_struct, ObjIndex, rhs_all, residual_flag)
%%%
%
% find the Lagrange multipliers through the factorization: corresponds to lexlse.ObjectiveSensitivity
% (for a single objective)
%

    if nargin < 4
	residual_flag = 0;
    end

    %% -------------------------------------
    %% input handling
    %% -------------------------------------
    obj_info = lexqr_struct.obj_info;
    lexqr    = lexqr_struct.lexqr;
    nVar     = lexqr_struct.nVar;
    hh_scalars = lexqr_struct.hh_scalars;
    %% -------------------------------------

    nLambda = 0; % Number of constraints that may influence objective ObjIndex
    nRank   = 0; % Total rank of objectives before objective ObjIndex
    for k=1:ObjIndex-1
	nLambda = nLambda + obj_info(k).dim;
	nRank   = nRank + obj_info(k).rank;
    end
    nLambda = nLambda + obj_info(ObjIndex).dim;

    Lambda = zeros(nLambda,1);
    rhs    = rhs_all(1:nRank);

    %% -------------------------------------

    FirstRowIndex = obj_info(ObjIndex).first_row_index;
    FirstColIndex = obj_info(ObjIndex).first_col_index;
    ObjDim        = obj_info(ObjIndex).dim;
    ObjRank       = obj_info(ObjIndex).rank;
    ColDim        = FirstColIndex-1; % ColDim will be used to indicate column dimension (I use it for clarity)

    if residual_flag
	Lambda = segment_assign(lexqr_struct.residual_mu{ObjIndex}, Lambda, FirstRowIndex, ObjDim);
    else
	%% Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed

	%% copy only what is needed to compute the residual v = A*x-b (i.e., -y_hat)
	tmp    = segment(lexqr(:,nVar+1), FirstRowIndex+ObjRank, ObjDim-ObjRank);
	Lambda = segment_assign(-tmp, Lambda, FirstRowIndex+ObjRank, ObjDim-ObjRank);

	%% compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
	tmp    = segment(Lambda, FirstRowIndex, ObjDim);
	hh_seq = householderSequence(block(lexqr, FirstRowIndex, FirstColIndex, ObjDim, ObjRank), ...
				     segment(hh_scalars, FirstRowIndex, ObjDim));
	tmp    = applyOnTheLeft(hh_seq,tmp);
	Lambda = segment_assign(tmp, Lambda, FirstRowIndex, ObjDim);
    end

    if (ObjIndex>1) % the first objective has only Lagrange multipliers equal to the optimal residual

        %% e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}

	tmp1 = block(lexqr, FirstRowIndex, 1, ObjDim, ColDim);
	tmp2 = segment(Lambda, FirstRowIndex, ObjDim);

	%%keyboard
	%%rhs = head_assign(-tmp1'*tmp2, rhs, ColDim);
	rhs = head_accumulate(-tmp1'*tmp2, rhs, ColDim);

	for k=ObjIndex-1:-1:1 % for all objectives before ObjIndex

            FirstRowIndex = obj_info(k).first_row_index;
	    FirstColIndex = obj_info(k).first_col_index;
	    ObjDim        = obj_info(k).dim;
	    ObjRank       = obj_info(k).rank;
	    ColDim        = FirstColIndex-1;

            %% Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed

	    tmp    = segment(rhs, FirstColIndex, ObjRank);
            Lambda = segment_assign(tmp, Lambda, FirstRowIndex, ObjRank);

            %% apply Q_k' on the left
	    tmp    = segment(Lambda, FirstRowIndex, ObjDim);
	    hh_seq = householderSequence(block(lexqr, FirstRowIndex, FirstColIndex, ObjDim, ObjRank), ...
					 segment(hh_scalars, FirstRowIndex, ObjDim));
	    tmp    = applyOnTheLeft(hh_seq,tmp);
	    Lambda = segment_assign(tmp, Lambda, FirstRowIndex, ObjDim);

	    tmp1 = block(lexqr, FirstRowIndex, 1, ObjDim, ColDim);
	    tmp2 = segment(Lambda, FirstRowIndex, ObjDim);

	    rhs = head_accumulate(-tmp1'*tmp2, rhs, ColDim);
	end

    end

%%%EOF
