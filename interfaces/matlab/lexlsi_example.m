% test 1
fprintf('=======================================\n     test 01\n')
clear lexobj options;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, w] = lexlsi(lexobj)

% test 2
fprintf('=======================================\n     test 02\n')
clear lexobj options;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, w] = lexlsi(lexobj,options)


% test 3
fprintf('=======================================\n     test 03\n')
clear lexobj options;
lexobj(1).A = rand(6,5);
lexobj(1).ub = rand(6,1);
lexobj(1).lb = lexobj(1).ub - rand(6,1);

lexobj(2).A = rand(4,5);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);
lexobj(2).c = [0; 1; 2; 0];


[x, info, w] = lexlsi(lexobj)



% test 4
fprintf('=======================================\n     test 04\n')
clear lexobj options;
options.enable_simple_bounds = 1;
lexobj(1).A = [1; 4; 5];
lexobj(1).ub = rand(3,1);
lexobj(1).lb = lexobj(1).ub - rand(3,1);
lexobj(1).c = [0; 1; 2];

lexobj(2).A = rand(4,6);
lexobj(2).ub = rand(4,1);
lexobj(2).lb = lexobj(2).ub - rand(4,1);

[x, info, w, active_set] = lexlsi(lexobj,options)
