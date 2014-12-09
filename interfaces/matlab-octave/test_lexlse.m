test_counter = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvarlist = ['clearvarlist', setdiff(who, {'test_counter'})];
clear(clearvarlist{:});

test_counter = test_counter + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options1.enable_fixed_variables = 1;
options1.regularization = [0.0, 0.01];
options1.regularizationType = 1;

options2.enable_fixed_variables = 1;
options2.regularization = [0.1, 0.01];
options2.regularizationType = 1;

options3.enable_fixed_variables = 0;
options3.regularization = [0.0, 0.01];
options3.regularizationType = 1;


num_var = 8;

lexobj12(1).A = [3;4;7]; % indices of the fixed variables
lexobj12(1).b = rand(3, 1);

lexobj12(2).A = rand(7, num_var);
lexobj12(2).b = rand(7, 1);


lexobj3(1).A = eye(num_var);
lexobj3(1).A = lexobj3(1).A(lexobj12(1).A, :);
lexobj3(1).b = lexobj12(1).b;

lexobj3(2).A = lexobj12(2).A;
lexobj3(2).b = lexobj12(2).b;


try
    [X1] = lexlse(lexobj12, options1);
    [X2] = lexlse(lexobj12, options2);
    [X3] = lexlse(lexobj3, options3);

    if (norm(X1 - X2) ~= 0) || (norm(X1 - X3) > 1e-10)
        fprintf(['test ', num2str(test_counter, '%04d'), ' ..... Failed\n']);
    else
        fprintf(['test ', num2str(test_counter, '%04d'), ' ..... OK\n']);
    end
catch
    fprintf(['test ', num2str(test_counter, '%04d'), ' ..... Failed\n']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
