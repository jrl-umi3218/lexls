%
% Copyright 2013-2021 INRIA
%

function [check_flag, out] = compare_results(obj,x1,x2,tol)
%
%
%

if nargin < 4
  tol = 1e-12;
end

nObj = length(obj);
check_flag = 0;

% -----------------------------------------------------------

for i=1:nObj
  if size(obj(1).A,2) == 1 && i == 1 % fixed variables
    r1{i} = x1(obj(i).A) - obj(i).b;
    r2{i} = x2(obj(i).A) - obj(i).b;
  else
    r1{i} = obj(i).A*x1 - obj(i).b;
    r2{i} = obj(i).A*x2 - obj(i).b;
  end

  err_residual(i) = norm(r1{i}-r2{i});

  if err_residual(i) > tol
    check_flag = 1;
  end
end

err_x = norm(x1-x2);
if err_x > tol
  check_flag = 1;
end

% -----------------------------------------------------------

out.r1 = r1;
out.r2 = r2;

out.err_residual = err_residual;
out.err_x        = err_x;

%%%EOF