function export_hierarchy(obj, file_name, active_set_guess, solution_guess, solution)
    HIERARCHY_NONE                  = 0;   %
    HIERARCHY_EQUALITIES            = 100; % equality constraints
    HIERARCHY_INEQUALITIES          = 200; % inequality constraints
    HIERARCHY_INEQUALITIES_WITH_AS  = 210; % inequality constraints with active set guess

    TYPES_OF_OBJECTIVES_SIMPLE      = 100; % simple constraints (fixed variables / simple bounds);
    TYPES_OF_OBJECTIVES_GENERAL     = 200; % general constraints

    % -------------------------------------------------------------------------

    % define identifiers
    NUMBER_OF_VARIABLES         = '#nVar';
    NUMBER_OF_OBJECTIVES        = '#nObj';
    NUMBER_OF_CONSTRAINTS       = '#nCtr';
    HIERARCHY_TYPE              = '#HierType';
    TYPES_OF_OBJECTIVES         = '#ObjType';
    OBJECTIVE_DATA              = '#OBJECTIVE';
    SOLUTION_GUESS              = '#SolGuess';
    SOLUTION                    = '#Solution';

    % -------------------------------------------------------------------------


    if (nargin < 5)
        solution = [];
    end

    if (nargin < 4)
        solution_guess = [];
    end

    if (nargin < 3)
        active_set_guess = [];
    end


    % -------------------------------------------------------------------------
    number_of_objectives = length(obj);
    number_of_variables = 0;
    number_of_constraints_per_level = [];

    hierarchy_type = HIERARCHY_NONE;
    types_of_objectives = [];
    % -------------------------------------------------------------------------



    % -------------------------------------------------------------------------
    if (isfield(obj(1), 'ub'))
        if (isempty(active_set_guess))
            hierarchy_type = HIERARCHY_INEQUALITIES;
        else
            hierarchy_type = HIERARCHY_INEQUALITIES_WITH_AS;
        end
    else
        hierarchy_type = HIERARCHY_EQUALITIES;
    end


    for i = 1:number_of_objectives;
        if (i == 1) && (size(obj(i).A, 2) == 1);
            types_of_objectives = [types_of_objectives, TYPES_OF_OBJECTIVES_SIMPLE];
        else
            if (i == 1) || ((i == 2) && (types_of_objectives(1) == TYPES_OF_OBJECTIVES_SIMPLE));
                number_of_variables = size(obj(i).A, 2);
            end

            types_of_objectives = [types_of_objectives, TYPES_OF_OBJECTIVES_GENERAL];
        end

        number_of_constraints_per_level = [number_of_constraints_per_level, size(obj(i).A, 1)];
    end
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    fileID = fopen(file_name,'w');


    fprintf(fileID, '# Exported at (%s) \n', datestr(clock));

    fprintf(fileID, ['\n', HIERARCHY_TYPE,'\n']);
    fprintf(fileID, '%d ', hierarchy_type);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NUMBER_OF_VARIABLES,'\n']);
    fprintf(fileID, '%d ', number_of_variables);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NUMBER_OF_OBJECTIVES,'\n']);
    fprintf(fileID, '%d ', number_of_objectives);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NUMBER_OF_CONSTRAINTS,'\n']);
    fprintf(fileID, '%d ', number_of_constraints_per_level);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', TYPES_OF_OBJECTIVES,'\n']);
    fprintf(fileID, '%d ', types_of_objectives);
    fprintf(fileID, '\n');


    if (hierarchy_type >= 100) && (hierarchy_type < 200)
        for i = 1:number_of_objectives
            fprintf(fileID,['\n', OBJECTIVE_DATA,' %d\n'],i-1);

            for j = 1:number_of_constraints_per_level(i);
                if (types_of_objectives(i) == TYPES_OF_OBJECTIVES_SIMPLE)
                    fprintf(fileID,'% d % 1.18e ', [obj(i).A(j), obj(i).b(j)]);
                else
                    fprintf(fileID,'% 1.18e ', [obj(i).A(j,:), obj(i).b(j)]);
                end
                fprintf(fileID,'\n');
            end
        end
    end

    if (hierarchy_type >= 200) && (hierarchy_type < 300)
        for i = 1:number_of_objectives
            fprintf(fileID, ['\n', OBJECTIVE_DATA,' %d\n'], i-1);

            for j = 1:number_of_constraints_per_level(i);
                if (types_of_objectives(i) == TYPES_OF_OBJECTIVES_SIMPLE)
                    fprintf(fileID,'% d % 1.18e % 1.18e ', [obj(i).A(j), obj(i).lb(j), obj(i).ub(j)]);
                else
                    fprintf(fileID,'% 1.18e ', [obj(i).A(j,:), obj(i).lb(j), obj(i).ub(j)]);
                end

                if (hierarchy_type == HIERARCHY_INEQUALITIES_WITH_AS)
                    fprintf(fileID,'% d ', active_set_guess{i}(j));
                end
                fprintf(fileID,'\n');
            end
        end
    end

    if (~isempty(solution_guess))
        fprintf(fileID, ['\n', SOLUTION_GUESS,'\n']);
        fprintf(fileID, '% 1.18e \n', solution_guess);
        fprintf(fileID, '\n');
    end

    if (~isempty(solution))
        fprintf(fileID, ['\n', SOLUTION,'\n']);
        fprintf(fileID, '% 1.18e \n', solution);
        fprintf(fileID, '\n');
    end

    fprintf(fileID,'\n');

    fclose(fileID);
end
