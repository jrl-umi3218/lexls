%
% Copyright 2013-2021 INRIA
%

function out = tail(vector,dim)
%%%
%
% returns the last dim entries of a vector
%

    if dim > length(vector)
	fprintf('dim > length(vector) \n')
	keyboard
    end

    out = vector(end-dim+1:end);

%%%EOF
