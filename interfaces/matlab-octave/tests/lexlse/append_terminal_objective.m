function [obj, options] = append_terminal_objective(obj, options)
%
%
%

if options.get_least_norm_solution
  nObj = length(obj);
  n    = size(obj(2).A,2);
  
  obj(nObj+1).A = eye(n);
  obj(nObj+1).b = zeros(n,1);
  
  options.regularization = [options.regularization, 0];
end

%%%EOF4