%
% Copyright 2013-2021 INRIA
%

function matrix = block_assign(B, matrix, row_start, col_start, row_dim, col_dim)
%%%
%
% assigns a block of a given matrix
%

    row_ind = row_start:row_start+row_dim-1;
    col_ind = col_start:col_start+col_dim-1;

    matrix(row_ind, col_ind) = B;

%%%EOF
