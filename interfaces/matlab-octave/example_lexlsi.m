% test 1
fprintf('=======================================\n     test 01\n')
clear;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, as, w] = lexlsi(lexobj)

% test 2
fprintf('=======================================\n     test 02\n')
clear;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, as, w] = lexlsi(lexobj, options)


% test 3
fprintf('=======================================\n     test 03\n')
clear;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set = {};
active_set{1} = [];
active_set{2} = [0; 1; 2; 0];


[x, info, as, w] = lexlsi(lexobj, [], active_set)



% test 4
fprintf('=======================================\n     test 04\n')
clear;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set_guess = {};
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

[x, info, active_set, w] = lexlsi(lexobj, options, active_set_guess)



% test 5
fprintf('=======================================\n     test 05\n')
clear;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set_guess = {};
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

x0 = zeros(6,1);

[x, info, active_set, w] = lexlsi(lexobj, options, active_set_guess, x0)



% test 6
fprintf('=======================================\n     test 06\n')
clear;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

x0 = zeros(6,1);

[x, info, active_set, w] = lexlsi(lexobj, options, [], x0)



% test 7
fprintf('=======================================\n     test 07\n')
clear;

options.enable_simple_bounds = 1;
options.regularization = [0.1, 0.1];
options.regularizationType = 1;

lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

x0 = zeros(6,1);

[x, info, active_set, w] = lexlsi(lexobj, options, [], x0)



% test 8
fprintf('=======================================\n     test 08\n')
clear;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set_guess = {};
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

x0 = zeros(6,1);

w = active_set_guess;

[x, info, active_set, w] = lexlsi(lexobj, options, active_set_guess, x0, w)



% test 9
fprintf('=======================================\n     test 09\n')
clear;
options.enable_simple_bounds = 1;
options.regularization = [0, 0.03];
options.regularizationType = 1;
options.variable_regularization_factor = 0.1;

lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set_guess = {};
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

x0 = zeros(6,1);

w = active_set_guess;

[x, info, active_set, w, debug] = lexlsi(lexobj, options, active_set_guess, x0, w)
