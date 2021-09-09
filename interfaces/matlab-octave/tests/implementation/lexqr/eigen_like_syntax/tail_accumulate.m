%
% Copyright 2013-2021 INRIA
%

function vector = tail_accumulate(v,vector,dim)
%%%
%
% accumulates to the tail of a given vector
%

    if dim > length(vector)
	fprintf('dim > length(vector) \n')
	keyboard
    end

    if dim ~= length(v)
	fprintf('dim ~= length(v) \n')
	keyboard
    end

    vector(end-dim+1:end) = vector(end-dim+1:end) + v;

%%%EOF
