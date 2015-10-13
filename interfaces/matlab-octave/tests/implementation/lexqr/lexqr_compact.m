function lexqr_struct = lexqr_compact(lexqr_struct)
%%%
%
% extract a compact lexqr
% (requires RT and LQ in lexqr_struct, generate using lexqr_explicit.m)
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
    R        = lexqr_struct.R;
    Q        = lexqr_struct.Q;
    %% -------------------------------------

    Qr = zeros(nCtr,nRank);
    Qs = zeros(nCtr,nCtr-nRank);
    Rr = zeros(nRank,nVar+1); %% note: Rs = 0 and is not stored

    rank_counter     = 1;
    singular_counter = 1;
    for ObjIndex=1:nObj

	FirstRowIndex = obj_info(ObjIndex).first_row_index;
	FirstColIndex = obj_info(ObjIndex).first_col_index;
	ObjDim          = obj_info(ObjIndex).dim;
	ObjRank         = obj_info(ObjIndex).rank;

	ind_rank_full        = FirstRowIndex:FirstRowIndex+ObjRank-1;
	ind_rank_compact     = rank_counter:rank_counter+ObjRank-1;
	ind_singular_full    = FirstRowIndex+ObjRank:FirstRowIndex+ObjDim-1;
	ind_singular_compact = singular_counter:singular_counter+ObjDim-ObjRank-1;

	Rr(ind_rank_compact,FirstColIndex:end) = R(ind_rank_full,FirstColIndex:end);
	Qr(:,ind_rank_compact)     = Q(:,ind_rank_full);
	Qs(:,ind_singular_compact) = Q(:,ind_singular_full);

	rank_counter     = rank_counter + ObjRank;
	singular_counter = singular_counter + ObjDim - ObjRank;
    end

    lexqr_struct.Rr = Rr;
    lexqr_struct.Qr = Qr; % full-rank part
    lexqr_struct.Qs = Qs; % singular part

%%%EOF
