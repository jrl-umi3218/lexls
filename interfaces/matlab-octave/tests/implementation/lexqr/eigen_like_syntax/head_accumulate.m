function vector = head_accumulate(v,vector,dim)
%%%
%
% accumulate to the head of a vector
%

    if dim > length(vector)
	fprintf('dim > length(vector) \n')
	keyboard
    end

    if dim ~= length(v)
	fprintf('dim ~= length(v) \n')
	keyboard
    end

    if dim > 0
	vector(1:dim) = vector(1:dim) + v;
    end

%%%EOF
