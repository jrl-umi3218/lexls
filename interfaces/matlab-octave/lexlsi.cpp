#include <vector>
#include <cmath>

#include <stdint.h>
#include <mex.h>
#include <typedefs.h>
#include <lexlsi.h>


#include "lexls_common.h"

const int MIN_INPUTS = 1;
const int MAX_INPUTS = 5;
const int MIN_OUTPUTS = 1;
const int MAX_OUTPUTS = 4;
const int MIN_NUMBER_OF_FIELDS_IN_OBJ = 3;


enum ConstraintActivationType
{
    BOUNDS_INACTIVE = 0,
    LOWER_BOUND_ACTIVE = 1,
    UPPER_BOUND_ACTIVE = 2,
    EQUALITY_CONSTRAINT = 3
};



mxArray * formInfoStructure (
        const LexLS::TerminationStatus lexlsi_status, 
        const int num_activations, 
        const int num_deactivations)
{
    mxArray * info_struct;

    int num_info_fields = 3; 
    const char *info_field_names[] = {
        "status",
        "number_of_activations",
        "number_of_deactivations",
    }; 


    info_struct = mxCreateStructMatrix(1, 1, num_info_fields, info_field_names);


    ReturnStatus status;
    switch (lexlsi_status)
    {
        case LexLS::PROBLEM_SOLVED:
            status = STATUS_OK;
            break;
        case LexLS::PROBLEM_SOLVED_CYCLING_HANDLING:
            status = STATUS_CYCLING_HANDLING;
            break;
        case LexLS::MAX_NUMBER_OF_ITERATIONS_EXCEDED:
            status = STATUS_MAXITER;
            break;
        default:
            status = STATUS_UNKNOWN_FAILURE;
            break;
    }
    mxArray *info_status = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    ((INT32_T *) mxGetData (info_status))[0] = status;
    mxSetField (info_struct, 0, "status", info_status);

    mxArray *info_activations = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    ((INT32_T *) mxGetData (info_activations))[0] = num_activations;
    mxSetField (info_struct, 0, "number_of_activations", info_activations);

    mxArray *info_deactivations = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    ((INT32_T *) mxGetData (info_deactivations))[0] = num_deactivations;
    mxSetField (info_struct, 0, "number_of_deactivations", info_deactivations);

    return (info_struct);
}


/**
 * @brief The main function.
 */
void mexFunction( int num_output, mxArray *output[],
                  int num_input, const mxArray *input[])
{
    checkInputOutput(num_output, output, num_input, input, 
                    MIN_OUTPUTS, MAX_OUTPUTS, MIN_INPUTS, MAX_INPUTS,
                    MIN_NUMBER_OF_FIELDS_IN_OBJ);

// parameters of the solver
    LexLS::LexLSIParameters lexlsi_parameters;

    lexlsi_parameters.setDefaults();


    bool    is_simple_bounds_handling_enabled = false;

    bool                is_regularization_set = false;
    std::vector<double> regularization;



// parse parameters

    if (num_input >= 2)
    {
        const mxArray *options_struct = input[1];

        if ((options_struct != NULL) && (!mxIsEmpty (options_struct)))
        {
            // ================================================
            // parse options

            getOptionDouble(&lexlsi_parameters.tolLinearDependence, 
                            options_struct, 
                            "tolLinearDependence");

            getOptionDouble(&lexlsi_parameters.tolWrongSignLambda, 
                            options_struct, 
                            "tolWrongSignLambda");

            getOptionDouble(&lexlsi_parameters.tolCorrectSignLambda, 
                            options_struct, 
                            "tolCorrectSignLambda");



            getOptionInteger(   &lexlsi_parameters.max_number_of_iterations, 
                                options_struct, 
                                "max_iterations");



            getOptionBool(  &lexlsi_parameters.CyclingHandling, 
                            options_struct, 
                            "cycling_handling");

            getOptionInteger(   &lexlsi_parameters.cycling_max_counter, 
                                options_struct, 
                                "cycling_max_counter");


            getOptionDouble(&lexlsi_parameters.cycling_relax_step, 
                            options_struct, 
                            "cycling_relax_step");



            is_regularization_set = getOptionArray( regularization, 
                                                    mxGetNumberOfElements (input[0]),
                                                    options_struct, 
                                                    "regularization");

            getOptionBool(  &is_simple_bounds_handling_enabled,
                            options_struct, 
                            "enable_simple_bounds");

            int regularization_type = 0;
            if (getOptionInteger( &regularization_type, 
                              options_struct, 
                              "regularizationType"))
            {
                lexlsi_parameters.regularizationType = static_cast <LexLS::RegularizationType> (regularization_type);
            }


            getOptionInteger(   &lexlsi_parameters.regularizationMaxIterCG, 
                                options_struct, 
                                "regularizationMaxIterCG");


            getOptionBool(  &lexlsi_parameters.realSensitivityResidual, 
                            options_struct, 
                            "realSensitivityResidual");


            // ================================================
            // check provided options

            /// @todo This check should probably go to the solver interface
            if ( 
                    ((lexlsi_parameters.regularizationType == LexLS::REGULARIZATION_NONE) && (is_regularization_set)) ||
                    ((lexlsi_parameters.regularizationType != LexLS::REGULARIZATION_NONE) && (!is_regularization_set))
               )
            {
                std::cout << lexlsi_parameters.regularizationType << std::endl;
                std::cout << is_regularization_set << std::endl;
                mexErrMsgTxt("Both regularization type and regularization factors must be specified.");
            }
        }
    }


// parse objectives
    const mxArray * objectives = input[0];

    LexLS::Index num_obj = mxGetNumberOfElements (objectives);

    std::vector<mxArray *> constraints;

    std::vector<LexLS::ObjectiveType> obj_type;
    std::vector<LexLS::Index> num_constr;
    std::vector<LexLS::Index> simple_bounds_indicies;


    obj_type.resize(num_obj);
    num_constr.resize(num_obj);
    constraints.resize(num_obj);

    int total_num_constr = 0;
    int first_normal_obj_index = 0;


    if (is_simple_bounds_handling_enabled)
    { // simple bounds

        mxArray *A = getObjectiveMatrix(objectives, 0, "A");
        mxArray *lb = getObjectiveMatrix(objectives, 0, "lb");
        mxArray *ub = getObjectiveMatrix(objectives, 0, "ub");


    // check A and b
        num_constr[0] = mxGetM(A);
        total_num_constr += num_constr[0];

        checkInputMatrixSize(A, num_constr[0], 1, 0, "A");
        checkInputMatrixSize(lb, num_constr[0], 1, 0, "lb");
        checkInputMatrixSize(ub, num_constr[0], 1, 0, "ub");


    // form[lb, ub]
    // remember objective
        constraints[0] = catenateMatrices(lb, ub);
        obj_type[0] = LexLS::SIMPLE_BOUNDS_OBJECTIVE_HP;

    // copy indices
        int num_simple_bounds = mxGetM(A);

        for (int i = 0; i < num_simple_bounds; ++i)
        {
            simple_bounds_indicies.push_back(static_cast <int> (round(mxGetPr(A)[i])) - 1);
        }

    // ok
        first_normal_obj_index = 1;

        failIfTrue(num_obj == 1, "Problems consisting of one level of simple bounds are not supported.");
    }


    LexLS::Index num_var = 0;

    for (int i = first_normal_obj_index; i < num_obj; ++i)
    {
        mxArray *A = getObjectiveMatrix(objectives, i, "A");
        mxArray *lb = getObjectiveMatrix(objectives, i, "lb");
        mxArray *ub = getObjectiveMatrix(objectives, i, "ub");

    // check A and b
        num_constr[i] = mxGetM(A);
        total_num_constr += num_constr[i];

        if (i == first_normal_obj_index)
        {
            num_var = mxGetN(A);
        }

        checkInputMatrixSize(A, num_constr[i], num_var, i, "A");
        checkInputMatrixSize(lb, num_constr[i], 1, i, "lb");
        checkInputMatrixSize(ub, num_constr[i], 1, i, "ub");

    // form Ab = [A, b]
    // remember objective
        constraints[i] = catenateMatrices(A, lb, ub);
        obj_type[i] = LexLS::DEFAULT_OBJECTIVE;
    }


// process active set guess
    std::vector<mxArray *> active_set;
    active_set.resize(num_obj);

    if (num_input >= 2)
    {
        const mxArray *active_set_cell = input[2];

        if ((active_set_cell != NULL) && (!mxIsEmpty (active_set_cell)))
        {
            failIfTrue (!mxIsCell(active_set_cell), "Active set must be of 'cell' type.");
            failIfTrue (mxGetNumberOfElements(active_set_cell) != num_obj, "Wrong dimention of the active set.");

            for (int i = 0; i < num_obj; ++i)
            {
                mxArray *c = mxGetCell (active_set_cell, i);
                if ((c == NULL) || (mxIsEmpty (c)))
                {
                    active_set[i] = NULL;
                }
                else
                {
                    failIfTrue (!mxIsDouble(c), "Active set flags must be of 'double' type.");
                    checkInputMatrixSize(c, num_constr[i], 1, "active set");
                    active_set[i] = c;
                }
            }
        }
    }


// process solution guess
    const mxArray *x0 = NULL;
    bool is_initial_guess_set = false;

    if (num_input >= 3)
    {
        x0 = input[3];

        if ((x0 != NULL) && (!mxIsEmpty (x0)) )
        {
            failIfTrue (!mxIsDouble(x0), "Initial guess must be of 'double' type.");

            checkInputMatrixSize(x0, num_var, 1, "x0");
            is_initial_guess_set = true;
        }
    }


// process residual guess
    std::vector<mxArray *> residuals;
    residuals.resize(num_obj);

    if (num_input >= 5)
    {
        const mxArray *residuals_cell = input[4];

        if ((residuals_cell != NULL) && (!mxIsEmpty (residuals_cell)))
        {
            failIfTrue (!mxIsCell(residuals_cell), "Residuals must be of 'cell' type.");
            failIfTrue (mxGetNumberOfElements(residuals_cell) != num_obj, "Wrong dimention of the residuals.");

            for (int i = 0; i < num_obj; ++i)
            {
                mxArray *c = mxGetCell (residuals_cell, i);
                if ((c == NULL) || (mxIsEmpty (c)))
                {
                    residuals[i] = NULL;
                }
                else
                {
                    failIfTrue (!mxIsDouble(c), "Residuals must be of 'double' type.");
                    checkInputMatrixSize(c, num_constr[i], 1, "residuals");
                    residuals[i] = c;
                }
            }
        }
    }


// instantiate LexLS
    LexLS::TerminationStatus status;
    LexLS::LexLSI lexlsi(num_var, num_obj, num_constr.data(), obj_type.data());


// initialize solver
    try
    {
        lexlsi.setParameters(lexlsi_parameters);


        // regularization
        if (is_regularization_set)
        {
            for (int i = 0; i < num_obj; ++i)
            {
                lexlsi.setRegularization(i, regularization[i]);
            }
        }

        // constraints
        if (is_simple_bounds_handling_enabled)
        {
            lexlsi.setData(
                    0, 
                    simple_bounds_indicies.data(), 
                    Eigen::Map<LexLS::MatrixType>(
                        mxGetPr(constraints[0]), 
                        num_constr[0], 
                        2));
        }
        for (int i = first_normal_obj_index; i < num_obj; ++i)
        {
            lexlsi.setData(
                    i, 
                    Eigen::Map<LexLS::MatrixType>(
                        mxGetPr(constraints[i]), 
                        num_constr[i], 
                        num_var + 2));
        }


        // set initial guess
        if (is_initial_guess_set)
        {
            /// @todo Map does not work for some reason
            //lexlsi.set_x0(Eigen::Map<LexLS::dVectorType>(mxGetPr(x0), num_var));
            LexLS::dVectorType x0_in(num_var);

            for (int i = 0; i < num_var; ++i)
            {
                x0_in(i) = (static_cast <const double *> (mxGetPr(x0)))[i];
            }

            lexlsi.set_x0(x0_in);
        }


        // set residuals
        for (int i = 0; i < num_obj; ++i)
        {
            if (residuals[i] == NULL)
            {
                // no residual
                continue;
            }
            else
            {
                LexLS::dVectorType residual(num_constr[i]);

                for (int j = 0; j < num_constr[i]; ++j)
                {
                    residual(j) = (static_cast <const double *> (mxGetPr(residuals[i])))[j];
                }

                lexlsi.set_w(i, residual);
            }
        }


        // activate constraints
        for (int i = 0; i < num_obj; ++i)
        {
            if (active_set[i] == NULL)
            {
                // nothing to activate;
                continue;
            }
            else
            {
                for (int j = 0; j < num_constr[i]; ++j)
                {
                    switch (static_cast <int> (round(mxGetPr(active_set[i])[j])))
                    {
                        case BOUNDS_INACTIVE:
                            // nothing to activate
                            break;
                        case LOWER_BOUND_ACTIVE:
                            lexlsi.api_activate(i, j, LexLS::LOWER_BOUND);
                            break;
                        case UPPER_BOUND_ACTIVE:
                            lexlsi.api_activate(i, j, LexLS::UPPER_BOUND);
                            break;
                        case EQUALITY_CONSTRAINT:
                            lexlsi.api_activate(i, j, LexLS::EQUALITY_CONSTRAINT);
                            break;
                        default:
                            mexWarnMsgTxt("Wrong format of 'c' matrix. Assuming inactive constraint.");
                            break;
                    }
                }
            }
        }

    // solve the problem
        status = lexlsi.solve();
    }
    catch (std::exception &e)
    {
        mexErrMsgTxt(e.what());
    }


// output solution
    try
    {
        LexLS::dVectorType& x = lexlsi.get_x();
        output[0] = mxCreateDoubleMatrix(num_var, 1, mxREAL);
        double *x_out = mxGetPr(output[0]);
        for (int i = 0; i < num_var; ++i)
        {
            x_out[i] = x(i);
        }
    }
    catch (std::exception &e)
    {
        mexErrMsgTxt(e.what());
    }

// output additional information
    if (num_output >= 2)
    {
        LexLS::Index num_activations;
        LexLS::Index num_deactivations;
       
        try
        {
            num_activations = lexlsi.getAddCount();
            num_deactivations = lexlsi.getRemoveCount();
        }
        catch (std::exception &e)
        {
            mexErrMsgTxt(e.what());
        }
        output[1] = formInfoStructure(status, num_activations, num_deactivations);
    }

// output residual
    if (num_output >= 3)
    {
        output[2] = mxCreateCellMatrix(num_obj, 1);
        for (int i = 0; i < num_obj; ++i)
        {

            try
            {
                LexLS::dVectorType& w = lexlsi.getResidual(i);
                mxArray * wi = mxCreateDoubleMatrix(num_constr[i], 1, mxREAL);
                for (int j = 0; j < num_constr[i]; ++j)
                {
                    mxGetPr(wi)[j] = w(j);
                }
            
                mxSetCell(output[2], i, wi);
            }
            catch (std::exception &e)
            {
                mexErrMsgTxt(e.what());
            }
        }
    }

// output active set
    if (num_output >= 4)
    {
        std::vector<LexLS::ConstraintType>  lexlsi_active_constraints;

        output[3] = mxCreateCellMatrix(num_obj, 1);

        for (int i = 0; i < num_obj; ++i)
        {
            try
            {
                lexlsi_active_constraints.clear();
                lexlsi.getActiveCtr(i, lexlsi_active_constraints);
            }
            catch (std::exception &e)
            {
                mexErrMsgTxt(e.what());
            }

            mxArray * active_constraints = mxCreateDoubleMatrix(lexlsi_active_constraints.size(), 1, mxREAL);
            for (unsigned int j = 0; j < lexlsi_active_constraints.size(); ++j)
            {
                switch (lexlsi_active_constraints[j])
                {
                    case LexLS::LOWER_BOUND:
                        mxGetPr(active_constraints)[j] = LOWER_BOUND_ACTIVE;
                        break;
                    case LexLS::UPPER_BOUND:
                        mxGetPr(active_constraints)[j] = UPPER_BOUND_ACTIVE;
                        break;
                    case LexLS::EQUALITY_CONSTRAINT:
                        mxGetPr(active_constraints)[j] = EQUALITY_CONSTRAINT;
                        break;
                    default:
                        mxGetPr(active_constraints)[j] = BOUNDS_INACTIVE;
                        break;
                }
            }
            mxSetCell(output[3], i, active_constraints);
        }
    }
}
