function vector = head_assign(v,vector,dim)
%%%
%
% assigns to the head of a vector
%

    if dim > length(vector)
	fprintf('dim > length(vector) \n')
	keyboard
    end

    if dim ~= length(v)
	fprintf('dim ~= length(v) \n')
	keyboard
    end

    vector(1:dim) = v;

%%%EOF
