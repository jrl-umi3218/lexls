%
% Copyright 2013-2021 INRIA
%

function obj = define_problem(n, m, r, fixed_variables)
%
% generate a random problem
%
%
% Input:
% ------
% n     - dimension of x
% m(i)  - number of constraints in level i
% r(i)  - rank that obj(i).A contributes to [obj(1).A; ...; obj(i-1).A]
% fixed_variables - define fixed variables at first level
%
% Output:
% -------
% obj   - problem definition
%

if nargin < 4
  fixed_variables = 0;
  if nargin < 3
    r = m;
  end
end

C = [];
for i=1:length(m)

  if fixed_variables && i == 1
    r(1) = m(1);

    obj(i).A = randperm(n,m(1))'; % generate m(1) unique integers in [1,n]
    obj(i).b = randn(m(i),1);

    ind      = obj(1).A;
    A        = zeros(m(1),n);
    A(:,ind) = eye(length(ind));

    C = [C;A];
  else
    obj(i).A = randn(m(i),sum(m(1:i-1)) + r(i)) * [C;randn(r(i),n)];
    obj(i).b = randn(m(i),1);

    C = [C;obj(i).A];
    %% scale (for many levels the norm of the matrices becomes large)
    s = max(max(abs(C)));
    if s > 1
	C = C/s;
    end
  end

end

%%%EOF
