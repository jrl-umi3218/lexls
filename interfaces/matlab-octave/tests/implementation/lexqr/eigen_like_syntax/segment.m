%
% Copyright 2013-2021 INRIA
%

function out = segment(vector, ind_start, dim)
%%%
%
% returns a segment of a vector
%

    out = vector(ind_start:ind_start+dim-1);

%%%EOF
