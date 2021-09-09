%
% Copyright 2013-2021 INRIA
%

function lexqr_struct = lexqr_explicit(lexqr_struct)
%%%
%
% form the orthonormal matrices explicitly
%

    %% -------------------------------------
    %% input handling
    %% -------------------------------------
    obj_info = lexqr_struct.obj_info;
    lexqr    = lexqr_struct.lexqr;
    nRank    = sum([lexqr_struct.obj_info(:).rank]);
    nCtr     = lexqr_struct.nCtr;
    nVar     = lexqr_struct.nVar;
    nObj     = lexqr_struct.nObj;
    hh_scalars = lexqr_struct.hh_scalars;
    %% -------------------------------------

    AccumulateObjDim = 1;
    Q = zeros(nCtr,nCtr);
    R = zeros(nRank,nVar+1);
    for ObjIndex=1:nObj
	FirstRowIndex = obj_info(ObjIndex).first_row_index;
	FirstColIndex = obj_info(ObjIndex).first_col_index;
	ObjDim        = obj_info(ObjIndex).dim;
	ObjRank       = obj_info(ObjIndex).rank;

	%% R
	ind_row_full = FirstRowIndex:FirstRowIndex+ObjDim-1;
	R(ind_row_full,FirstColIndex:end) = triu(lexqr(ind_row_full,FirstColIndex:end));

	hh_seq = householderSequence(block(lexqr, FirstRowIndex, FirstColIndex, ObjDim, ObjRank), ...
				     segment(hh_scalars, FirstRowIndex, ObjDim));

	%% Q
	Qk = applyOnTheLeft(hh_seq, eye(ObjDim)); % form Q explicitly
	Q  = block_assign(Qk,Q,FirstRowIndex,FirstRowIndex,ObjDim,ObjDim);

	%% left blocks
	RemainingRows = nCtr-(FirstRowIndex+ObjDim-1);
	Lk = block(lexqr, FirstRowIndex+ObjDim, FirstColIndex, RemainingRows, ObjRank);
	Q  = block_assign(Lk, Q, FirstRowIndex+ObjDim, AccumulateObjDim, RemainingRows, ObjRank);

	AccumulateObjDim = AccumulateObjDim + ObjDim;
    end

    lexqr_struct.Q = Q;
    lexqr_struct.R = R;

%%%EOF
