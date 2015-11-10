function [x,L,v,obj,d] = wset_solve(wset, obj, options)
%%%
%
%
%

    nObj = length(obj);

    for i=1:nObj
	nCtr = length(wset{i});
	nActive(i) = 0;

	J = [];
	k = 1;
	for j=1:nCtr

	    if wset{i}(j) == 1               % lower bound
		obj(i).ub(j) = obj(i).lb(j);

		nActive(i) = nActive(i) + 1;
	    elseif wset{i}(j) == 2           % upper bound
		obj(i).lb(j) = obj(i).ub(j);

		nActive(i) = nActive(i) + 1;
	    elseif wset{i}(j) == 3           % equality
		%% do nothing
	    else                             % inactive
		J(k) = j;
		k = k+1;
	    end
	end
	obj(i).A(J,:) = [];
	obj(i).lb(J)  = [];
	obj(i).ub(J)  = [];
    end

    [x, info, v, as, d] = lexlsi(obj,options);

    L = [];
    for i=1:nObj
	L = [L;d.lambda{i}];
    end

%%%EOF
