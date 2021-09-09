%
% Copyright 2013-2021 INRIA
%

function [x,L,v,obj_eq,d,A,x_eq] = wset_solve_order(active_ctr, obj, options)
%%%
%
%
%

    nObj = length(obj);
    nVar = size(obj(1).A,2);
    nActiveCtr = length(active_ctr);

    for i=1:nObj
	obj_eq(i).A  = zeros(0,nVar);
	obj_eq(i).lb = zeros(0,1);
	obj_eq(i).ub = zeros(0,1);
	obj_eq(i).b = zeros(0,1);
    end

    A = [];
    b = [];

    counter = ones(1,nObj);
    for i=1:nActiveCtr

	obj_index = active_ctr(i).obj_index+1;
	ctr_index = active_ctr(i).ctr_index+1;
	ctr_type  = active_ctr(i).ctr_type;
	k         = counter(obj_index);

	A = [A;obj(obj_index).A(ctr_index,:)];

	if ctr_type == 1     % lower bound

	    obj_eq(obj_index).A(k,:)  = obj(obj_index).A(ctr_index,:);
	    obj_eq(obj_index).lb(k,:) = obj(obj_index).lb(ctr_index);
	    obj_eq(obj_index).ub(k,:) = obj(obj_index).lb(ctr_index);

	    obj_eq(obj_index).b(k,:) = obj(obj_index).lb(ctr_index);

	    b = [b;obj(obj_index).lb(ctr_index)];

	elseif ctr_type == 2 % upper bound

	    obj_eq(obj_index).A(k,:)  = obj(obj_index).A(ctr_index,:);
	    obj_eq(obj_index).lb(k,:) = obj(obj_index).ub(ctr_index);
	    obj_eq(obj_index).ub(k,:) = obj(obj_index).ub(ctr_index);

	    obj_eq(obj_index).b(k,:) = obj(obj_index).ub(ctr_index);

	    b = [b;obj(obj_index).ub(ctr_index)];

	elseif ctr_type == 3 % equality (use upper bound by convention)

	    obj_eq(obj_index).A(k,:)  = obj(obj_index).A(ctr_index,:);
	    obj_eq(obj_index).lb(k,:) = obj(obj_index).ub(ctr_index);
	    obj_eq(obj_index).ub(k,:) = obj(obj_index).ub(ctr_index);

	    obj_eq(obj_index).b(k,:) = obj(obj_index).ub(ctr_index);

	    b = [b;obj(obj_index).ub(ctr_index)];
	else

	    disp('ERROR: we shoudl not be here')

	end

	counter(obj_index) = counter(obj_index) + 1;
    end

    A = [A,b];

    [x, info, v, as, d] = lexlsi(obj_eq,options);
    x_eq = lexlsi(obj_eq,options);

    L = [];
    for i=1:nObj
	L = [L;d.lambda{i}];
    end

%%%EOF
