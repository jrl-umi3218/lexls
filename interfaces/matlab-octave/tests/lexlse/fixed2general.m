%
% Copyright 2013-2021 INRIA
%

function [obj, options] = fixed2general(obj, options)
%
% from obj(1) with fixed variables generate obj(1) with general equality constraints
%

if nargin < 2
  options = [];
end

[m1,n1] = size(obj(1).A);
if n1 == 1
  ind      = obj(1).A;
  n        = size(obj(2).A,2);
  A        = zeros(m1,n);
  A(:,ind) = eye(length(ind));
  obj(1).A = A;
  end

  options.enable_fixed_variables = 0;
end

%%%EOF
