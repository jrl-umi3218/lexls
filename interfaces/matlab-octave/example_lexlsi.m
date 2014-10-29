% test 1
fprintf('=======================================\n     test 01\n')
clear;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, w] = lexlsi(lexobj)

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

[x, info, w] = lexlsi(lexobj, [], [], options)


% test 3
fprintf('=======================================\n     test 03\n')
clear;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

active_set = cell();
active_set{1} = [];
active_set{2} = [0; 1; 2; 0];


[x, info, w] = lexlsi(lexobj, active_set)



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

active_set_guess = cell();
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

[x, info, w, active_set] = lexlsi(lexobj, active_set_guess, [], options)



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

active_set_guess = cell();
active_set_guess{1} = [0; 1; 2];
active_set_guess{2} = [];

x0 = zeros(6,1);

[x, info, w, active_set] = lexlsi(lexobj, active_set_guess, x0, options)



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

[x, info, w, active_set] = lexlsi(lexobj, [], x0, options)



% test 7
fprintf('=======================================\n     test 07\n')
clear;

options.enable_simple_bounds = 1;
options.regularization = [0.1, 0.1];

lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

x0 = zeros(6,1);

[x, info, w, active_set] = lexlsi(lexobj, [], x0, options)