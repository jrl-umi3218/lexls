% test 1
fprintf('=======================================\n     test 01\n')
clear lexobj options;
lexobj(1).A = rand(3,5);
lexobj(1).b = rand(3,1);

lexobj(2).A = rand(4,5);
lexobj(2).b = rand(4,1);

[x, info, w] = lexlse(lexobj)


% test 2 (pseudoinverse)
fprintf('=======================================\n     test 02\n')
clear lexobj options;
lexobj(1).A = rand(3,5);
lexobj(1).b = rand(3,1);

lexobj(2).A = eye(5);
lexobj(2).b = zeros(5,1);

[x, info, w] = lexlse(lexobj)

x1 = pinv(lexobj(1).A) * lexobj(1).b;
err = norm(x - x1)


% test 3 (fixed variables)
fprintf('=======================================\n     test 03\n')
clear lexobj options;
options.enable_fixed_variables = 1;
lexobj(1).A = [3;4;7]; % indices of the fixed variables
lexobj(1).b = rand(3, 1);

lexobj(2).A = rand(7, 8);
lexobj(2).b = rand(7, 1);

[x, info, w] = lexlse(lexobj, options)
lexobj(1)


% test 4 (tolerance)
fprintf('=======================================\n     test 04\n')
clear lexobj options;
options.linear_dependence_tolerance = 10^-3;
lexobj(1).A = rand(5, 8);
lexobj(1).b = rand(5, 1);

lexobj(2).A = rand(6, 8);
lexobj(2).b = rand(6, 1);

[x, info, w] = lexlse(lexobj, options)


% test 5 (info)
fprintf('=======================================\n     test 05\n')
clear lexobj options;
lexobj(1).A = rand(5, 8);
lexobj(1).b = rand(5, 1);

lexobj(2).A = rand(6, 8);
lexobj(2).b = rand(6, 1);

[x, info, w] = lexlse(lexobj)
