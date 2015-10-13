%%%
%
% verify implementation
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')

clear;clc

%% -----------------------------------------------------------------------

n = 15;
m = 5*ones(4,1);
r = m-2;

e = [];
for i=1:100

    %% ===================================================================

    lexqr_struct = define_problem(n,m,r);

    %% ===================================================================

    [x,info,v,as,d] = lexlsi(lexqr_struct.obj, lexqr_struct.options);

    %% ===================================================================

    lexqr_struct = lexqr_form(lexqr_struct);
    lexqr_struct = lexqr_solve(lexqr_struct);
    lexqr_struct = lexqr_residual(lexqr_struct);
    lexqr_struct = lexqr_explicit(lexqr_struct);
    lexqr_struct = lexqr_compact(lexqr_struct);
    lexqr_struct = lexqr_lambda(lexqr_struct);
    vec          = lexqr_residual_compact(lexqr_struct);
    L            = lexqr_lambda_compact(lexqr_struct);

    %% ===================================================================

    e_name{1} = 'factor';
    e(i,1) = norm(lexqr_struct.lexqr - d.lexqr);

    e_name{2} = 'x';
    e(i,2) = norm(x-lexqr_struct.x);

    e_name{3} = 'residual';
    e(i,3) = norm(vec-lexqr_struct.v);

    L0 = [];
    for k=1:lexqr_struct.nObj
	L0 = [L0;d.lambda{k}];
    end
    e_name{4} = 'lambda';
    e(i,4) = norm(L-lexqr_struct.lambda);

    fprintf('i = %4d (',i)
    fprintf(' %e ', e(i,:))
    fprintf(')\n')

    out_probelm{i} = lexqr_struct.obj;

    %% ===================================================================

end

%% -----------------------------------------------------------------------

fprintf('\n')

for i=1:size(e,2)
    [min_val(i), min_ind(i)] = min(e(:,i));
    [max_val(i), max_ind(i)] = max(e(:,i));
    mean_val(i)              = mean(e(:,i));
end

fprintf(2,'---------------------------------------------------------------------------------\n')
fprintf('       ')
for i=1:length(e_name)
    fprintf('%d: [%s]     ', i, e_name{i})
end
fprintf('\n')
fprintf(2,'---------------------------------------------------------------------------------\n')
fprintf(2,' min: ')
fprintf(' %e ', min_val)
fprintf('\n')

fprintf(2,' max: ')
fprintf(' %e ', max_val)
fprintf('\n')
fprintf(2,'---------------------------------------------------------------------------------\n')

%%%EOF
