function vector = segment_assign(v,vector,ind_start,dim)
%%%
%
% assigns a vector to a segment of a vector
%

    vector(ind_start:ind_start+dim-1) = v;

%%%EOF
