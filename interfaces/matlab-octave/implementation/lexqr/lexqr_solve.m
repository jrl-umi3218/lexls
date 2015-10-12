function lexqr_struct = lexqr_solve(lexqr_struct)
%%%
%
% find x given the lexqr decomposition: corresponds to lexlse.solve()
%

    %% -------------------------------------
    %% input handling
    %% -------------------------------------
    obj_info = lexqr_struct.obj_info;
    lexqr    = lexqr_struct.lexqr;
    P        = lexqr_struct.P;
    nVar     = lexqr_struct.nVar;
    nObj     = lexqr_struct.nObj;
    %% -------------------------------------

    x = zeros(nVar,1);
    AccumulatedRanks = 0;

    for k=nObj:-1:1

        ObjRank = obj_info(k).rank;
        if (ObjRank > 0)

            x_k = segment(lexqr(:,nVar+1), obj_info(k).first_row_index, ObjRank);

	    if (AccumulatedRanks > 0) % Do not enter here the first time ObjRank != 0

		tmp = block(lexqr, ...
			    obj_info(k).first_row_index, obj_info(k+1).first_col_index, ...
			    ObjRank, AccumulatedRanks);

		x_k = x_k - tmp * segment(x, obj_info(k+1).first_col_index, AccumulatedRanks);
	    end

		tmp = block(lexqr, ...
			    obj_info(k).first_row_index, obj_info(k).first_col_index, ...
			    ObjRank, ObjRank);

		x_k = triu(tmp)\x_k;
		x = segment_assign(x_k, x, obj_info(k).first_col_index, ObjRank);

		AccumulatedRanks = AccumulatedRanks + ObjRank;
	end
    end

    x = P*x;

    lexqr_struct.x = x;

%%%EOF
