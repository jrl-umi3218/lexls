function A = swap(A,c1,c2)
%%%
%
% swap A(:,c1) and A(:,c2)
%

    n = size(A,2);
    if (c1 > n) || (c2 > n)
	fprintf('(c1 > n) || (c2 > n) \n')
	keyboard
    end

    tmp = A(:,c2);
    A(:,c2) = A(:,c1);
    A(:,c1) = tmp;

%%%EOF
