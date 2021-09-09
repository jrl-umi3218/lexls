%
% Copyright 2013-2021 INRIA
%

function [A,tau,beta,h] = makeHouseholderInPlace(A,row_start,row_dim,col)
%%%
%
% Given a vector: segment(A(:,col),row_start,row_dim),
% store the essential part of the corresponding Householder vector in:
% segment(A(:,col),row_start+1,row_dim-1)
%
% The entire Householder vector is given by h = [1;essential] and can
% be used to form the Householder reflection H = eye(row_dim)-tau*h*h'
%
% Note: This corresponds exactly to the implementation in Eiegn (for real data)
%
% -------------------------------------------------------------------
% Example:
% -------------------------------------------------------------------
% A = randn(6,4);
%
% col       = 3;
% row_start = 2;
% row_dim   = 5;
%
% [A,tau,beta,h] = makeHouseholderInPlace(A,row_start,row_dim,col);
% H = eye(row_dim)-tau*h*h';
%
% H*A(2:end,:)
% -------------------------------------------------------------------

    v = segment(A(:,col), row_start,row_dim);

    essential = v(2:end);
    vec_tail  = essential;

    if length(v) == 1
	tailSqNorm = 0;
    else
	tailSqNorm = squaredNorm(vec_tail);
    end
    c0 = v(1);

    if(tailSqNorm == 0)
	tau  = 0;
	beta = c0;
	essential(:) = 0;
    else
	 beta = sqrt(c0*c0 + tailSqNorm);
	 if (c0 >= 0)
	     beta = -beta;
	 end
	 essential = vec_tail / (c0 - beta);
	 tau       = (beta - c0) / beta;
    end

    A = block_assign(essential, A, row_start+1, col, row_dim-1, 1);
    h = [1;essential]; % output for convenience

%%%EOF
