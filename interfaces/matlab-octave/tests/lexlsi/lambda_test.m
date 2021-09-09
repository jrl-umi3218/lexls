%
% Copyright 2013-2021 INRIA
%

%% ---------------------------------------------------------------
%% a test example from A. Escande
%% ---------------------------------------------------------------
%% (x_1 = 1) > (2*x_2 = 1) > ... > (n*x_n = 1) > (sum(x_i) = 1)
%%
%% lambda = zeros(n+1,n+1);
%% lambda(:,n+1) = [-w; -w/2; ...; -w/n; w];
%%
%% with w = sum(1./[2:n]);
%% ---------------------------------------------------------------

addpath('/Users/drdv/Work/Git/bip/soft/lexls/interfaces/matlab-octave')

clear;clc

n = 5;
w = sum(1./[2:n]);

lambda = zeros(n+1,n+1);
lambda(end,end) = w;
for k = 1:n
    obj(k).A  = zeros(1,n); obj(k).A(k) = k;
    obj(k).lb = 1;
    obj(k).ub = obj(k).lb;

    lambda(k,end) = -w/k;
end

k = n+1;
obj(k).A  = ones(1,n);
obj(k).lb = 1;
obj(k).ub = obj(k).lb;

%% ---------------------------------------------------------------

[x,info,as,v,d] = lexlsi(obj);

L = [];
for k = 1:n+1
    L = [L; d.lambda{k}];
end

norm(lambda - L)
