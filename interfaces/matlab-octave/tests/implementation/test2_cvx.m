%%%
%
%
%

addpath('./lexqr')
addpath('./utility')
addpath('./lexlse_dual')

clear;clc

nObj = 3;
nVar = 20;
m    = 5*ones(nObj,1);
r    = m-2;
mu   = 0.5*ones(nObj,1);

nCtr = sum(m);

lexqr_struct = define_problem(nVar,m,r);
obj = lexqr_struct.obj;

[Ab, m] = obj2array(obj);
A = Ab(:,1:end-1);

%% ---------------------------------------------------------
%% CVX (testing with 3 objectives)
%% ---------------------------------------------------------

cvx_precision best
cvx_clear

cvx_begin quiet
  variables x1(nVar) r1(m(1))
  dual variables L11

  minimize(0.5*r1'*r1 + 0.5*mu(1)^2*x1'*x1)
  subject to
    L11: r1 == obj(1).A*x1 - obj(1).b
cvx_end

fprintf('primal(1): cvx_status = %s\n',cvx_status);

cvx_begin quiet
  variables x2(nVar) r2(m(2))
  dual variables L21 L22

  minimize(0.5*r2'*r2 + 0.5*mu(2)^2*x2'*x2)
  subject to
    L21: r1 == obj(1).A*x2 - obj(1).b
    L22: r2 == obj(2).A*x2 - obj(2).b
cvx_end

fprintf('primal(2): cvx_status = %s\n',cvx_status);

cvx_begin quiet
  variables x3(nVar) r3(m(3))
  dual variables L31 L32 L33

  minimize(0.5*r3'*r3 + 0.5*mu(3)^2*x3'*x3)
  subject to
    L31: r1 == obj(1).A*x3 - obj(1).b
    L32: r2 == obj(2).A*x3 - obj(2).b
    L33: r3 == obj(3).A*x3 - obj(3).b
cvx_end

fprintf('primal(3): cvx_status = %s\n',cvx_status);

LL = zeros(nCtr,nObj);
i=1; LL(1:sum(m(1:i)),i) = L11;
i=2; LL(1:sum(m(1:i)),i) = [L21;L22];
i=3; LL(1:sum(m(1:i)),i) = [L31;L32;L33];

A'*LL + [mu(1)^2*x1,mu(2)^2*x2,mu(3)^2*x3]

[~, Ld, ~] = lexlse_dual(obj,mu);
disp(' LL - Ld (could be different):')
disp(LL - Ld)

%%%EOF
