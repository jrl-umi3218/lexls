%
% Copyright 2013-2021 INRIA
%

function x = seq_lexls(obj_in,mu,flag_basic)
%
% sequence of equality constrainted problems
%

nObj = length(obj_in);

for i=1:nObj
  A{i} = obj_in(i).A;
  b{i} = obj_in(i).b;
end

nVar = size(A{1},2);

if (length(mu) ~= nObj)
    mu = ones(1,nObj)*mu;
end

if nargin < 3
    flag_basic = 0;
end

% ------------------------------------------

I = eye(nVar);
z = zeros(nVar,1);

% ------------------------------------------

obj(1).A = [A{1};mu(1)*I];
obj(1).b = [b{1};z];

x = lexlse(obj);

if flag_basic
    if nObj == 1 % in order to get a basic solution at the end
        obj(1).A = [A{1}];
        obj(1).b = [A{1}*x];

        x = lexlse(obj);
    end
end

for k = 2:nObj

  C = []; c = [];
  for i = 1:k-1
    C = [C;A{i}];
    c = [c;A{i}*x];
  end

  obj(1).A = C;
  obj(1).b = c;
  obj(2).A = [A{k};mu(k)*I];
  obj(2).b = [b{k};z];

  x = lexlse(obj);

  if flag_basic
      if k == nObj % in order to get a basic solution at the end
          obj(1).A = C;
          obj(1).b = c;
          obj(2).A = [A{k}];
          obj(2).b = [A{k}*x];

          x = lexlse(obj);
      end
  end

end

%%%EOF
