function matrix = applyOnTheLeft(hh_seq, matrix)
%%%
%
% apply a sequence of Householder transformations on the left of "matrix"
%

    [n,m] = size(matrix);
    for i=length(hh_seq):-1:1
	matrix = applyHouseholderOnTheLeft(hh_seq(i).v, hh_seq(i).s, matrix, i, 1, n-i+1, m);
    end

%%%EOF
