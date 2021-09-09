%
% Copyright 2013-2021 INRIA
%

function [x, exitflag, info] = qpsolver(solver, x0, H, h, lb, ub, A, b, Ain, Alb, Aub)
    if (isempty(solver))
        solver = 'qpoases'
    end

    h = h(:);
    b = b(:);

    if (nargin == 8)
        Ain = [];
        Alb = [];
        Aub = [];
    else
        Alb = Alb(:);
        Aub = Aub(:);
    end


    t0 = clock ();
    switch solver
    case 'qpoases'
        options = qpOASES_options('default');
        options.enableRamping = 0;
        %options.enableRegularisation = 1;
        %options.printLevel = 3;
        options.initialStatusBounds = 0;

        options.enableRegularisation = 0;
        options.numRegularisationSteps = 0;

        options.enableFlippingBounds = 0;
        options.enableNZCTests = 0;

        options.enableEqualities = 1;

        %options.terminationTolerance = 2.22044604925031e-07;
        %options.terminationTolerance = 2.22044604925031e-13;
        %options.boundTolerance = 2.22044604925031e-13;

        % manual regularization seems to be more reliable
        H = H + eye(size(H))*options.epsRegularisation;

        A = [A; Ain];
        Alb = [b; Alb];
        Aub = [b; Aub];

	%% drdv: put to suppress warnings
	options.hessianType        = [];
	options.x0                 = [];
	options.guessedWorkingSetB = [];
	options.guessedWorkingSetC = [];
	options.R                  = [];

        [x,fval,exitflag,iter,lambda] = qpOASES(H, h, A, lb, ub, Alb, Aub, x0, options);
        info.fval = fval;
        info.iter = iter;
        info.lambda = lambda;

    case 'qpoases-mpc'
        options = qpOASES_options( 'MPC' );
        options.enableDriftCorrection = 1;  % without this the error in constraint satisfaction is higher (^-12).
        %options.printLevel = 0;
        %options.printLevel = 3;


        A = [A; Ain];
        Alb = [b; Alb];
        Aub = [b; Aub];

        [x, fval, exitflag, iter, lambda] = qpOASES(H, h, A, lb, ub, Alb, Aub, x0, options);
        info.fval = fval;
        info.iter = iter;
        info.lambda = lambda;

    case 'qld'
        Ain = [-Ain; Ain];
        bin = [-Alb; Aub];
        [x, obj, info, lambda] = qld(x0, H, h, A, b, lb, ub, Ain, bin);
        info.lambda = lambda;
        exitflag = info.info;

    case 'qld-hightol'
        Ain = [-Ain; Ain];
        bin = [-Alb; Aub];
        options.tolerance = 1e-8;
        [x, obj, info, lambda] = qld(x0, H, h, A, b, lb, ub, Ain, bin, options);
        info.lambda = lambda;
        exitflag = info.info;


    case 'octave'
        if (isempty(Ain))
            [x, obj, info, lambda] = qp (x0, H, h, A, b, lb, ub);
        else
            [x, obj, info, lambda] = qp (x0, H, h, A, b, lb, ub, Alb, Ain, Aub);
        end
        info.lambda = lambda;
        exitflag = info.info;


    otherwise
        error('Unknown solver');
    end

    info.exec_time = etime (clock (), t0);
end
