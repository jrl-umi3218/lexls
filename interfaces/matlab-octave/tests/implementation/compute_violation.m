function [v] = compute_violation(lexobj, X)
    v = {};

    for i = 1:numel(lexobj)
        if (i == 1) && (size(lexobj(i).A, 2) == 1)
            % assume simple bounds
            Ax = X(lexobj(i).A);
        else
            Ax = lexobj(i).A*X;
        end

        ind_lb = Ax < lexobj(i).lb;
        ind_ub = Ax > lexobj(i).ub;

        v_i = zeros(size(lexobj(i).A,1), 1);

        v_i(ind_lb) = Ax(ind_lb) - lexobj(i).lb(ind_lb);
        v_i(ind_ub) = Ax(ind_ub) - lexobj(i).ub(ind_ub);

        v = [v, {v_i}];
    end
end
