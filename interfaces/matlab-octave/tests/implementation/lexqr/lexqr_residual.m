function lexqr_struct = lexqr_residual(lexqr_struct)
  %%%
  %
  % find the residual through the factorization: corresponds to lexlse.get_v()
  %

  %% -------------------------------------
  %% input handling
  %% -------------------------------------
  obj_info = lexqr_struct.obj_info;
  lexqr    = lexqr_struct.lexqr;
  nCtr     = lexqr_struct.nCtr;
  nVar     = lexqr_struct.nVar;
  nObj     = lexqr_struct.nObj;
  hh_scalars = lexqr_struct.hh_scalars;
  %% -------------------------------------

  v = zeros(nCtr,1);

  for ObjIndex=1:nObj

    FirstRowIndex = obj_info(ObjIndex).first_row_index;
    FirstColIndex = obj_info(ObjIndex).first_col_index;
    ObjDim        = obj_info(ObjIndex).dim;
    ObjRank       = obj_info(ObjIndex).rank;

    tmp = segment(lexqr(:,nVar+1), FirstRowIndex+ObjRank, ObjDim-ObjRank);
    v   = segment_assign(-tmp, v, FirstRowIndex+ObjRank, ObjDim-ObjRank);

    tmp    = segment(v, FirstRowIndex, ObjDim);
    hh_seq = householderSequence(block(lexqr, FirstRowIndex, FirstColIndex, ObjDim, ObjRank), ...
      segment(hh_scalars, FirstRowIndex, ObjDim));
    tmp    = applyOnTheLeft(hh_seq,tmp);
    v      = segment_assign(tmp, v, FirstRowIndex, ObjDim);
  end

  lexqr_struct.v = v;

  for ObjIndex=1:nObj

    FirstRowIndex = obj_info(ObjIndex).first_row_index;
    ObjDim        = obj_info(ObjIndex).dim;

    lexqr_struct.obj_info(ObjIndex).v_obj = segment(v,FirstRowIndex,ObjDim);
  end

  %%%EOF
