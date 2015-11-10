function [x,info_out,w] = lex_sequence(lexobj,options,lambda)
%
% with damping
%

% ------------------------------------------

nObj = length(lexobj);
nVar = size(lexobj(1).A,2);

if nargin < 3
  lambda = zeros(1,nObj);
end

for k=1:nObj
  A{k} = lexobj(k).A;
  lb{k} = lexobj(k).lb;
  ub{k} = lexobj(k).ub;
end

I = eye(nVar);
z = zeros(nVar,1);

% ------------------------------------------

obj(1).A  = [A{1};lambda(1)*I];
obj(1).lb = [lb{1};z];
obj(1).ub = [ub{1};z];

[x,info,w] = lexlsi(obj,options);
W{1} = w{1}(1:length(lb{1}));

status(1) = info.status;
nA(1) = info.number_of_activations;
nD(1) = info.number_of_deactivations;
nF(1) = info.number_of_factorizations;
nC(1) = info.cycling_counter;

for k = 2:nObj

  for i=1:k-1
    obj(i).A  = A{i};
    obj(i).lb = lb{i}+W{i};
    obj(i).ub = ub{i}+W{i};
  end

  obj(k).A  = [A{k};lambda(k)*I];
  obj(k).lb = [lb{k};z];
  obj(k).ub = [ub{k};z];

  [x,info,w] = lexlsi(obj,options);
  W{k} = w{k}(1:length(lb{k}));

  status(k) = info.status;
  nA(k) = info.number_of_activations;
  nD(k) = info.number_of_deactivations;
  nF(k) = info.number_of_factorizations;
  nC(k) = info.cycling_counter;
end

w = W;
info_out.status = status;
info_out.nA = nA;
info_out.nD = nD;
info_out.nF = nF;
info_out.nC = nC;

%%%EOF
