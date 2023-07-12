%
% Copyright 2013-2021 INRIA
%

function test_lexlse(m,n,r,options,tol,counter)
%
%
%

fprintf('\n--------------------------------------------\n');
fprintf('test %d: ', counter);

obj = define_problem(n, m, r, options.enable_fixed_variables);
x1 = lexlse(obj,options);

[obj, options_] = fixed2general(obj, options);
[obj, options_] = append_terminal_objective(obj, options_);
x2 = lexlse(obj, options_);

[check_flag, out] = compare_results(obj,x1,x2,tol);

if ~check_flag
  fprintf('OK \n');
  fprintf('--------------------------------------------\n');

  if isfield(options,'get_least_norm_solution')
    fprintf('get_least_norm_solution = %d \n', options.get_least_norm_solution);
  end

  if isfield(options,'enable_fixed_variables')
    fprintf('enable_fixed_variables  = %d \n', options.enable_fixed_variables);
  end

  if isfield(options,'regularization_type')
    fprintf('regularization_type     = %d \n', options.regularization_type);
  end

  if isfield(options,'regularization_factors')
    fprintf('regularization_factors  = [');
    fprintf(' %d ', options.regularization_factors);
    fprintf('] \n');
  end
else
  fprintf('FAIL \n');
  fprintf('--------------------------------------------\n');

  fprintf('err_x        = %e \n', out.err_x);
  fprintf('err_residual = [');
  fprintf(' %e ', out.err_residual);
  fprintf('] \n');

  keyboard
end

%%%EOF
