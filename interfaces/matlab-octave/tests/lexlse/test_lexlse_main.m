%
% Tikhonov regularization using qr(A)
%

addpath('../../')

format short
clear;clc

% ---------------------------------------------------------

n = 30;
m = [9,8,10,6];
r = [7,6,8,5];

tol = 1e-10;

if 1 % check vs. sequence

  load_labels
  
  % ---------------------------------------------------------
  options.get_least_norm_solution = 0;
  options.enable_fixed_variables  = 0;
  options.regularization_type     = REGULARIZATION_TIKHONOV;
  options.regularization_factors  = [1,2,3,4];
  % ---------------------------------------------------------
  
  obj = define_problem(n, m, r, options.enable_fixed_variables);
  x1 = lexlse(obj,options);
  
  [obj, options_] = fixed2general(obj, options);
  [obj, options_] = append_terminal_objective(obj, options_);
  x2 = seq_lexls(obj,options.regularization_factors,1);
  
  [check_flag, out] = compare_results(obj,x1,x2,tol);
  
  % the solution could be different 
  % (the reason is related to column permutation)
  %[x1,x2,x1-x2] 

  if norm(out.err_residual) > tol
    keyboard
  end

end

options2test = test_lexlse_define();

for i=1:length(options2test)
  test_lexlse(m,n,r,options2test(i),tol,i);
end

%%%EOF