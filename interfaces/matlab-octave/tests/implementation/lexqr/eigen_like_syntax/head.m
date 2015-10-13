function out = head(vector,dim)
%%%
%
% returns the first dim entries of a vector
%

    if dim > length(vector)
	fprintf('dim > length(vector) \n')
	keyboard
    end

    out = vector(1:dim);

%%%EOF
