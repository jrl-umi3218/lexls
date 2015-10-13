function out = householderSequence(V, s)
%%%
%
%
%

    out = [];
    for i = 1:size(V,2)
	out(i).v = V(i+1:end,i);
	out(i).s = s(i);
    end

%%%EOF
