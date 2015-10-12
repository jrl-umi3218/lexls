function lexqr_struct = lexqr_form(lexqr_struct)
%%%
%
% form lexqr: corresponds to lexlse.factorize()
%

    [lexqr, m] = obj2array(lexqr_struct.obj);

    %% todo: need a function to handle options in general
    if isfield(lexqr_struct,'options')
	tol_linear_dependence = lexqr_struct.options.tol_linear_dependence;
    else
	tol_linear_dependence = 1e-12;
    end

    %% --------------------------------------------------------
    %% initializations
    %% --------------------------------------------------------

    nVar = size(lexqr,2)-1; % the last column is assumed to be the rhs vector
    nObj = length(m);
    nCtr = sum(m);
    hh_scalars = zeros(1,nCtr);

    P = eye(nVar); % permutation matrix

    obj_info(1).first_row_index = 1;
    obj_info(1).dim             = m(1);
    obj_info(1).rank            = 0;
    for i=2:nObj
	 obj_info(i).first_row_index = obj_info(i-1).first_row_index + m(i-1);
	 obj_info(i).dim             = m(i);
	 obj_info(i).rank            = 0;
    end

    %% --------------------------------------------------------
    %% form lexqr
    %% --------------------------------------------------------

    RemainingColumns = nVar;
    ColIndex         = 1;
    for ObjIndex=1:nObj  % loop over all objectives

	FirstRowIndex                      = obj_info(ObjIndex).first_row_index;
	FirstColIndex                      = ColIndex;
	obj_info(ObjIndex).first_col_index = ColIndex;
        ObjDim                             = obj_info(ObjIndex).dim;

        for k=ColIndex:nVar % initially compute the norms of the columns
	    ColNorms(k) = squaredNorm(segment(lexqr(:,k),FirstRowIndex,ObjDim));
	end

        for counter=0:ObjDim-1 %  loop over all constraints in a given objective

            RowIndex      = FirstRowIndex + counter; % current row to process
	    RemainingRows = ObjDim        - counter; % remaining rows in the current objective

            [maxColNormValue, maxColNormIndex] = maxCoeff(tail(ColNorms,RemainingColumns));
	    maxColNormIndex = maxColNormIndex + ColIndex - 1;

            maxColNormValue = squaredNorm(segment(lexqr(:,maxColNormIndex),RowIndex,RemainingRows));
	    ColNorms(maxColNormIndex) = maxColNormValue;

            if (maxColNormValue < tol_linear_dependence)
		break;
	    end

            %% --------------------------------------------------------------------------
	    %% apply column permutation
	    %% --------------------------------------------------------------------------
	    column_permutations(ColIndex) = maxColNormIndex;
	    if (ColIndex ~= maxColNormIndex)
		lexqr    = swap(lexqr,ColIndex,maxColNormIndex);
		ColNorms = swap(ColNorms,ColIndex,maxColNormIndex);
	    end

            %% --------------------------------------------------------------------------
	    %% apply Householder transformations (on the RHS as well)
	    %% --------------------------------------------------------------------------
	    if (RemainingRows > 1)

		[lexqr,tau,PivotValue]   = makeHouseholderInPlace(lexqr,RowIndex,RemainingRows,ColIndex);
		lexqr(RowIndex,ColIndex) = PivotValue;

		essential = segment(lexqr(:,ColIndex), RowIndex+1,RemainingRows-1);

		%% apply transformation on the RHS as well
		lexqr = applyHouseholderOnTheLeft(essential, tau, lexqr, ...
						  RowIndex, ColIndex+1, RemainingRows, RemainingColumns);
        	hh_scalars(FirstRowIndex+counter) = tau;
	    end
	    %% --------------------------------------------------------------------------

            ColIndex = ColIndex + 1;
	    RemainingColumns = nVar - ColIndex + 1;

	    %% terminate the QR factorization (after the elimination step below, the lexqr is terminated as well)
	    if (RemainingColumns == 0)
		break;
	    end

            %% update our table of squared norms of the columns
	    %% (note that above ColIndex is incremented and RemainingColumns is decreased)
	    if (RemainingRows > 0)
		s = segment(lexqr(RowIndex,:), ColIndex,RemainingColumns);
		ColNorms = tail_accumulate(-s.*s,ColNorms,RemainingColumns);
	    end
	end

	%% Note that here ColIndex is the index of the next available variable
	ObjRank = ColIndex - FirstColIndex;
	obj_info(ObjIndex).rank = ObjRank; % store the rank

        %% -----------------------------------------------------------------------
	%% Gauss transformation
	%% -----------------------------------------------------------------------

        if (ObjIndex < nObj) % if there are objectives with lower priority
	    if (ObjRank > 0)   % if there are variables to eliminate

		FirstRowIndexNextObjective = FirstRowIndex + ObjDim;
		RemainingRows = nCtr-FirstRowIndexNextObjective + 1; % VARIABLE REDEFINITION: remaining rows after

                LeftBlock     = block(lexqr, ...
				      FirstRowIndexNextObjective, FirstColIndex, ...
				      RemainingRows, ObjRank);

		UpBlock       = block(lexqr, ...
				      FirstRowIndex, ColIndex, ...
				      ObjRank, RemainingColumns+1);

		TrailingBlock = block(lexqr, ...
				      FirstRowIndexNextObjective, ColIndex, ...
				      RemainingRows, RemainingColumns+1);

		TriangBlock   = block(lexqr, ...
				      FirstRowIndex, FirstColIndex, ...
				      ObjRank, ObjRank);

		LeftBlock     = (triu(TriangBlock)'\LeftBlock')';
		TrailingBlock = TrailingBlock - LeftBlock * UpBlock;

		lexqr = block_assign(LeftBlock, lexqr, ...
				     FirstRowIndexNextObjective, FirstColIndex, ...
				     RemainingRows, ObjRank);

		lexqr = block_assign(TrailingBlock, lexqr, ...
				     FirstRowIndexNextObjective,ColIndex, ...
				     RemainingRows,RemainingColumns+1);
	    end
	end

	if (RemainingColumns == 0) % ColIndex = nVar
	    %% Initialize some remaining fields of obj_info before we terminate (used later in ComputeLambda())
	    for k=ObjIndex+1:nObj
		%% of course, obj_info(>ObjIndex).rank = 0
		obj_info(k).first_col_index = obj_info(k-1).first_col_index + obj_info(k-1).rank;
	    end

	    break; % terminate the lexqr
	end
    end

    TotalRank = 0;
    for ObjIndex=1:nObj
	TotalRank = TotalRank + obj_info(ObjIndex).rank;
    end

    %% form the permutation matrix
    for k=1:TotalRank
	P = swap(P,k,column_permutations(k));
    end

    %% --------------------------------------------------------
    %% form output
    %% --------------------------------------------------------
    k = 1;
    vars_used = [];
    for i=1:nObj
	first_row = obj_info(i).first_row_index;
	first_col = obj_info(i).first_col_index;
	for j=1:obj_info(i).rank
	    %% -----------------------------------------------------------
	    obj_info(i).pivots(j) = lexqr(first_row+j-1,first_col+j-1);
	    %% -----------------------------------------------------------
	    var_i = find(P(:,k)); % variables dedicated to i-th level;
	    obj_info(i).variables(j) = var_i;

	    vars_used = [vars_used, var_i];
	    k = k + 1;
	    %% -----------------------------------------------------------
	end
    end

    vars_unused = setdiff(1:nVar,vars_used);

    lexqr_struct.lexqr      = lexqr;       % factorization
    lexqr_struct.P          = P;           % permutation matrix
    lexqr_struct.nVar       = nVar;        % number of variables
    lexqr_struct.nObj       = nObj;        % number of objective
    lexqr_struct.nCtr       = nCtr;        % total number of constraints
    lexqr_struct.dof        = vars_unused; % not used variables
    lexqr_struct.obj_info   = obj_info;    % info structure
    lexqr_struct.hh_scalars = hh_scalars;  % Householder scalars
    %% --------------------------------------------------------

%%%EOF
