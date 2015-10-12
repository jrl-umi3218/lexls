function matrix = applyHouseholderOnTheLeft(essential, tau, matrix, row_start, col_start, row_dim, col_dim)
%%%
%
% Follows Eigen's implementation
%

    A      = block(matrix, row_start, col_start, row_dim  , col_dim);
    bottom = block(     A,         2,         1, row_dim-1, col_dim);

    tmp = essential(:)' * bottom;
    tmp = tmp + A(1,:);
    A(1,:) = A(1,:) - tau * tmp;
    bottom = bottom - tau * essential * tmp;

    A      = block_assign(bottom,      A,         2,         1, row_dim-1, col_dim);
    matrix = block_assign(     A, matrix, row_start, col_start, row_dim  , col_dim);

%%%EOF
