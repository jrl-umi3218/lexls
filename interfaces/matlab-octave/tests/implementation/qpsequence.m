% Time-stamp: <2013-06-12 22:17:58 drdv>
% Sequence to QPs for solving a lexicographic LS problem.
% The k-th QP has the form
%
% minimize_{x,w_k}  norm(w_k,2)^2
% subject to        lb_k <= A_k*x - w_k <= ub_k
%                   + constraints from objectives with higher priority
function [x,w] = qpsequence(Obj, regularization, reglevel)
%    solver = 'qld-hightol';
    solver = 'qpoases';

    % -----------------------------

    if (nargin < 2)
        regularization = 0;
    end

    if (nargin < 3)
        reglevel = 1;
    end

    num_var = size(Obj(1).A,2);
    num_obj = length(Obj);
    for k = 1:num_obj
        Obj(k).nCtr = size(Obj(k).A,1);
    end

    for k = 1:num_obj
        A = []; lbA = []; ubA = [];
        for j = 1:k-1
            A   = [A; Obj(j).A, zeros(Obj(j).nCtr, Obj(k).nCtr)];
            lbA = [lbA; Obj(j).lb + w{j}];
            ubA = [ubA; Obj(j).ub + w{j}];
        end
        A   = [A;Obj(k).A, -eye(Obj(k).nCtr)];
        lbA = [lbA; Obj(k).lb];
        ubA = [ubA; Obj(k).ub];

        if (k < reglevel)
            regflag = 0;
        else
            regflag = 1;
        end

        H  = diag([regflag*regularization*ones(num_var,1); ones(Obj(k).nCtr,1)]);
        h  = zeros(num_var+Obj(k).nCtr,1);


        [x, exitflag] = qpsolver(solver, [], H, h, [], [], [], [], A, lbA, ubA);
        if (exitflag ~= 0)
            fprintf('qp: exitflag = %d\n',exitflag);
            keyboard
        end

        w{k} = x(num_var+1:num_var+Obj(k).nCtr);
    end

    x = x(1:num_var);
end
