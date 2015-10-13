function [x_star, Lambda, X] = lexlse_dual(obj,mu)
%
%
%

    %% ----------------------------------------------------------
    %% input handling
    %% ----------------------------------------------------------

    nVar = size(obj(1).A,2);
    nObj = length(obj);

    for i=1:nObj
	A{i} = obj(i).A;
	y{i} = obj(i).b;
	m(i) = size(obj(i).A,1);
    end

    %% ----------------------------------------------------------

    x_star = zeros(nVar,1);
    for i=1:nObj
	B =[];
	for i=1:i
	    B = [B, A{i}'];
	end
	B = [B;zeros(m(i),sum(m(1:i-1))), mu(i)*eye(m(i))];
	b = -mu(i)*[mu(i)*x_star; y{i} - A{i}*x_star];

	S{i}   = B;
	L{i}   = pinv(B)*b;
	x_star = -1/(mu(i)^2)*B(1:nVar,:)*L{i};
	X(:,i) = x_star;
	R{i}   = B*L{i} - b;
    end

    nCtr = sum(m);
    Lambda = zeros(nCtr,nObj);
    for i=1:nObj
	Lambda(1:sum(m(1:i)),i) = L{i};
    end

%%%EOF
