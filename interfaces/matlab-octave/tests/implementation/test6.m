%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%


addpath('/Users/drdv/Work/Installed/cpp/qpOASES-3.1.0/interfaces/matlab')
addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')
addpath('./lexlse_dual')
addpath('./wset')

clear;%clc
format long

%% ===================================================================

if 1
    load test_debug.mat

    a = [];
    for i=1:length(debug)
	a = [a;debug(i).info.status];
    end

    obj = debug(113).tasks;

    A = [];
    b = [];
    for i=1:length(obj)
	A = [A;obj(i).A];
	b = [b;obj(i).lb,obj(i).ub];
    end
    A = [A,b];

    for i=1:length(obj)
	ind{i} = [];
	k = 1;
	for j=1:size(obj(i).A,1)
	    if norm([obj(i).A(j,:)]) < 1e-12
		ind{i}(k) = j;
		k = k+1;
	    end
	end
    end

    %% remove zero constraints
    if 0
	for i=1:length(obj)
	    obj(i).A(ind{i},:) = [];
	    obj(i).lb(ind{i})  = [];
	    obj(i).ub(ind{i})  = [];
	end
    end

    %%return
else

    if 0
	n = 8;
	m = 5*ones(3,1);

	for i=1:length(m)
	    obj(i).A  = randn(m(i),n);
	    obj(i).lb = randn(m(i),1);
	    obj(i).ub = obj(i).lb + rand(m(i),1);
	end

    else

	load mitko.mat

	obj = test;

    end

    if 0
	x = lexlsi(obj);
	export_hierarchy(obj, '/Users/drdv/Desktop/lex/cpp/test_tmp.dat', [], [], x);
	return
    end


end

nObj = length(obj);
nVar = size(obj(1).A,2);

%% ===================================================================

options.max_number_of_factorizations = 113;

%{
options.tol_linear_dependence   = 1e-18;
options.tol_wrong_sign_lambda   = 1e-14;
options.tol_correct_sign_lambda = 1e-13;
options.tol_feasibility         = 1e-13;
%}

%{
options.cycling_handling_enabled = 1;
options.cycling_max_counter      = options.max_number_of_factorizations;
options.cycling_relax_step       = 1e-05;
%}

options.deactivate_first_wrong_sign = 0;

%% ===================================================================

%{
for i=1:length(obj)
    for j=1:size(obj(i).A,1)
	if norm([obj(i).A(j,:)]) < 1e-08
	    disp('AHA')
	    keyboard
	end
    end
end
%}

[x,info,v0,as,d] = lexlsi(obj, options);

disp('-----------------------------')
fprintf('status                   = %d \n', info.status)
fprintf('number_of_factorizations = %d \n', info.number_of_factorizations)
disp('-----------------------------')

export_hierarchy(obj, '/Users/drdv/Desktop/lex/cpp/test_113.dat');

nLog = length(d.working_set_log);
wset = wset_get(as, d.working_set_log, nLog);

%%[x1,L1,v1,obj1,d1] = wset_solve(wset, obj, options);
[x2,L2,v2,obj2,d2,A2,x_eq] = wset_solve_order(d.active_ctr, obj, options);

x0 = d.xStar;

a = size(d.active_ctr);
norm(d.data(1:a,:) - d2.data(1:a,:))
norm(d.lexqr(1:a,:) - d2.lexqr(1:a,:))
norm(x2 - d.xStar)

return

%{
export_hierarchy(obj, '/Users/drdv/Desktop/lex/cpp/test_113.dat', [], [], d.xStar);

nLog = length(d.working_set_log);
wset = wset_get(as, d.working_set_log, nLog);

[x1,L1,v1,obj1,d1] = wset_solve(wset, obj, options);
[x2,L2,v2,obj2,d2,A2,x_eq] = wset_solve_order(d.active_ctr, obj, options);

export_hierarchy(obj2, '/Users/drdv/Desktop/lex/cpp/test_eq.dat', [], [], x2);

x0 = d.xStar;

a = size(d.active_ctr);
norm(d.data(1:a,:) - d2.data(1:a,:))
norm(d.lexqr(1:a,:) - d2.lexqr(1:a,:))
%}
%%norm(d.data(1:a,:) - A2)

%%[x0,x1,x2,x0-x1,x0-x2,x1-x2]
%%fprintf('norm(x1-x2) = %e \n', norm(x1-x2));

%%[x0,x2,x0-x2]
%fprintf('norm(x0-x2) = %e \n', norm(x0-x2));

k = 1;
ctr = [];
for j=1:size(A2,1)
    if norm(A2(j,:)) == 0
	ctr(k) = j;
	k = k+1;
    end
end

length([ind{:}])
length(ctr)
size(A2,1) - length(ctr)

return

for i=1:nObj
    ev01(i) = norm(v0{i}) - norm(v1{i});
    ev02(i) = norm(v0{i}) - norm(v2{i});
    ev12(i) = norm(v1{i}) - norm(v2{i});
end
fprintf('ev01 = %e \n', norm(ev01));
fprintf('ev02 = %e \n', norm(ev02));
fprintf('ev12 = %e \n', norm(ev12));

%% ===================================================================

for i=1:nObj
    e_as(i) = norm(wset{i} - as{i});
end
fprintf('active_set_error = %e \n', norm(e_as));

%%%EOF
