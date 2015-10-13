function [Ab, m] = obj2array(obj)
%%%
%
%
%

    A = [];
    b = [];
    for i=1:length(obj)
	A = [A;obj(i).A];
	b = [b;obj(i).b];
	m(i) = length(obj(i).b);
    end
    Ab = [A,b];

%%%EOF
