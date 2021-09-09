%
% Copyright 2013-2021 INRIA
%

function out = block(matrix, row_start, col_start, row_dim, col_dim)
%%%
%
% returns a block of a given matrix
%

    row_ind = row_start:row_start+row_dim-1;
    col_ind = col_start:col_start+col_dim-1;

    out = matrix(row_ind, col_ind);

%%%EOF
