%%%
%
% Compute Lambda for regularized problem using the compact lexqr
%

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')
addpath(genpath('lexqr'))
addpath('./utility')

clear;clc

%% -----------------------------------------------------------------------

n = 15;
m = 5*ones(4,1);
r = m-2;

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

%% ===================================================================

b = randn(16,1);

RT = lexqr_struct.Rr(:,1:end-1)';
QT = lexqr_struct.Qr';
P  = lexqr_struct.P;
c  = inv(RT'*RT)*RT'*P'*lexqr_struct.x;

%Qr*l(k) = mu(k)*c;

%%%EOF
