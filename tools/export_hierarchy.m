function export_hierarchy(obj, file_name)
    HIERARCHY_NONE = 0;            % 
    HIERARCHY_EQUALITIES = 1;      % equality constraints
    HIERARCHY_INEQUALITIES = 2;    % inequality constraints

    OBJTYPE_SIMPLE  = 100;    % simple constraints (fixed variables / simple bounds);
    OBJTYPE_GENERAL = 200;    % general constraints

    % -------------------------------------------------------------------------

    % define identifiers
    NVAR = '#nVar';
    NOBJ = '#nObj';
    NCTR = '#nCtr';
    TYPE = '#Type';
    OBJTYPE = '#ObjType';
    DATA = '#OBJECTIVE';
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    number_of_levels = length(obj);
    number_of_variables = 0;
    number_of_constraints_per_level = [];

    hierarchy_type = HIERARCHY_NONE;
    types_of_objectives = [];
    % -------------------------------------------------------------------------



    % -------------------------------------------------------------------------
    if (isfield(obj(1), 'ub'))
        hierarchy_type = HIERARCHY_INEQUALITIES;
    else
        hierarchy_type = HIERARCHY_EQUALITIES;
    end


    for i = 1:number_of_levels;
        if (i == 1) && (size(obj(i).A, 2) == 1);
            types_of_objectives = [types_of_objectives, OBJTYPE_SIMPLE];
        else
            if (i == 1) || ((i == 2) && (types_of_objectives(1) == OBJTYPE_SIMPLE));
                number_of_variables = size(obj(i).A, 2);
            end

            types_of_objectives = [types_of_objectives, OBJTYPE_GENERAL];
        end

        number_of_constraints_per_level = [number_of_constraints_per_level, size(obj(i).A, 1)];
    end
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    fileID = fopen(file_name,'w');


    fprintf(fileID, '# Exported at (%s) \n', datestr(clock));

    fprintf(fileID, ['\n', TYPE,'\n']);
    fprintf(fileID, '%d ', hierarchy_type);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NVAR,'\n']);
    fprintf(fileID, '%d ', number_of_variables);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NOBJ,'\n']);
    fprintf(fileID, '%d ', number_of_levels);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', NCTR,'\n']);
    fprintf(fileID, '%d ', number_of_constraints_per_level);
    fprintf(fileID, '\n');

    fprintf(fileID, ['\n', OBJTYPE,'\n']);
    fprintf(fileID, '%d ', types_of_objectives);
    fprintf(fileID, '\n');


    if (hierarchy_type == HIERARCHY_EQUALITIES)
        for i = 1:number_of_levels
            fprintf(fileID,['\n', DATA,' %d\n'],i-1);

            for j = 1:number_of_constraints_per_level(i);
                fprintf(fileID,'% 1.15f ', [obj(i).A(j,:), obj(i).b]);
                fprintf(fileID,'\n');
            end
        end
    end

    if (hierarchy_type == HIERARCHY_INEQUALITIES)
        for i = 1:number_of_levels
            fprintf(fileID, ['\n', DATA,' %d\n'], i-1);

            for j = 1:number_of_constraints_per_level(i);
                fprintf(fileID,'% 1.15f ', [obj(i).A(j,:), obj(i).lb(j), obj(i).ub(j)]);
                fprintf(fileID,'\n');
            end
        end
    end

    fprintf(fileID,'\n');

    fclose(fileID);
end
