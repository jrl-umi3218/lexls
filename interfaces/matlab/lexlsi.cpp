#include <vector>
#include <cmath>

#include <mex.h>
#include <typedefs.h>
#include <lexlsi.h>


#include "lexls_common.h"

const int MIN_OUTPUTS = 1;
const int MAX_OUTPUTS = 4;
const int MIN_NUMBER_OF_FIELDS_IN_OBJ = 3;


enum ConstraintActivationType
{
    BOUNDS_INACTIVE = 0,
    LOWER_BOUND_ACTIVE = 1,
    UPPER_BOUND_ACTIVE = 2,
    EQUALITY_CONSTRAINT = 3,
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
    checkInputOutput(num_output, output, num_input, input, MIN_OUTPUTS, MAX_OUTPUTS, MIN_NUMBER_OF_FIELDS_IN_OBJ);

// parameters of the solver
    double linear_dependence_tolerance = 0.0;
    double deactivation_tolerance = 0.0;
    bool tolerances_are_given = false;

    bool simple_bounds_enabled = false;
    bool max_iterations_set = false;
    bool cycling_option_enabled = false;

    int max_iterations;
    int max_stall_iterations;

    if (num_input == 2)
    {
        mxArray * linear_dependence_tolerance_opt   = mxGetField (input[1], 0, "linear_dependence_tolerance");
        mxArray * deactivation_tolerance_opt        = mxGetField (input[1], 0, "deactivation_tolerance");

        if ((linear_dependence_tolerance_opt == NULL) || (deactivation_tolerance_opt == NULL))
        {
            //mexWarnMsgTxt("Some tolerances are not defined, the default ones are used.");
        }
        else
        {
            failIfTrue (!mxIsDouble(linear_dependence_tolerance_opt), "Tolerance must be of 'double' type.");
            failIfTrue (!mxIsDouble(deactivation_tolerance_opt     ), "Tolerance must be of 'double' type.");

            linear_dependence_tolerance     = mxGetPr(linear_dependence_tolerance_opt)[0];
            deactivation_tolerance          = mxGetPr(deactivation_tolerance_opt     )[0];

            tolerances_are_given = true;
        }


        mxArray * enable_simple_bounds_option = mxGetField (input[1], 0, "enable_simple_bounds");

        if (enable_simple_bounds_option != NULL) 
        {
            failIfTrue (!mxIsDouble(enable_simple_bounds_option), "Flag options must be of 'double' type.");
            if (*mxGetPr(enable_simple_bounds_option) != 0)
            {
                simple_bounds_enabled = true;
            }
        }

        mxArray * max_iterations_option = mxGetField (input[1], 0, "max_iterations");

        if (max_iterations_option != NULL) 
        {
            failIfTrue (!mxIsDouble(max_iterations_option), "Flag options must be of 'double' type.");
            max_iterations = static_cast <int> (round(*mxGetPr(max_iterations_option)));
            max_iterations_set = true;
        }

        mxArray * enable_cycling_handling_option = mxGetField (input[1], 0, "cycling_handling");

        if (enable_cycling_handling_option != NULL) 
        {
            failIfTrue (!mxIsDouble(enable_cycling_handling_option), "Flag options must be of 'double' type.");
            if (*mxGetPr(enable_cycling_handling_option) != 0)
            {
                max_stall_iterations = static_cast <int> (round(*mxGetPr(enable_cycling_handling_option)));
                cycling_option_enabled = true;
            }
        }

    }

// parse parameters
    LexLS::Index num_obj = mxGetNumberOfElements (input[0]);

    std::vector<mxArray *> constraints;

    std::vector<LexLS::ObjectiveType> obj_type;
    std::vector<LexLS::Index> num_constr;
    std::vector<LexLS::Index> simple_bounds_indicies;

    std::vector<mxArray *> active_set;

    obj_type.resize(num_obj);
    num_constr.resize(num_obj);
    active_set.resize(num_obj);
    constraints.resize(num_obj);

    int total_num_constr = 0;
    int first_normal_obj_index = 0;


    if (simple_bounds_enabled)
    { // simple bounds
        mxArray *A = mxGetField (input[0], 0, "A");
        mxArray *ub = mxGetField (input[0], 0, "ub");
        mxArray *lb = mxGetField (input[0], 0, "lb");
        mxArray *c = mxGetField (input[0], 0, "c");

    // check A and b
        checkInputMatrix(A, "A");
        checkInputMatrix(lb, "lb");
        checkInputMatrix(ub, "ub");

        num_constr[0] = mxGetM(A);
        total_num_constr += num_constr[0];

        checkInputMatrixSize(A, num_constr[0], 1, "A");
        checkInputMatrixSize(lb, num_constr[0], 1, "lb");
        checkInputMatrixSize(ub, num_constr[0], 1, "ub");


    // update active set
        if ((c == NULL) || (mxIsEmpty (c)))
        {
            active_set[0] = NULL;
        }
        else
        {
            failIfTrue (!mxIsDouble(c), "Matrix 'c' must be of 'double' type.");
            checkInputMatrixSize(c, num_constr[0], 1, "c");
            active_set[0] = c;
        }
          
    // form[lb, ub]
        mxArray *cat_input[3];
        mxArray *cat_output[1];

        mxArray *dimension = mxCreateDoubleMatrix (1, 1, mxREAL);
        *mxGetPr(dimension) = 2;

        cat_input[0] = dimension;
        cat_input[1] = lb;
        cat_input[2] = ub;

        if(mexCallMATLAB (1, cat_output, 3, cat_input, "cat") != 0)
        {
            mexErrMsgTxt("catenation of A and b failed!");
        }
        mxDestroyArray(dimension);

    // remember objective
        constraints[0] = cat_output[0];
        obj_type[0] = LexLS::SIMPLE_BOUNDS_OBJECTIVE_HP;

    // copy indices
        int num_simple_bounds = mxGetM(A);

        for (int i = 0; i < num_simple_bounds; ++i)
        {
            simple_bounds_indicies.push_back(static_cast <int> (round(mxGetPr(A)[i])) - 1);
        }

    // ok
        first_normal_obj_index = 1;
    }


    LexLS::Index num_var = 0;

    for (int i = first_normal_obj_index; i < num_obj; ++i)
    {
        mxArray *A = mxGetField (input[0], i, "A");
        mxArray *ub = mxGetField (input[0], i, "ub");
        mxArray *lb = mxGetField (input[0], i, "lb");
        mxArray *c = mxGetField (input[0], i, "c");

    // check A and b
        checkInputMatrix(A, "A");
        checkInputMatrix(lb, "lb");
        checkInputMatrix(ub, "ub");

        num_constr[i] = mxGetM(A);
        total_num_constr += num_constr[i];

        if (i == first_normal_obj_index)
        {
            num_var = mxGetN(A);
        }

        checkInputMatrixSize(A, num_constr[i], num_var, "A");
        checkInputMatrixSize(lb, num_constr[i], 1, "lb");
        checkInputMatrixSize(ub, num_constr[i], 1, "ub");

    // update active set
        if ((c == NULL) || (mxIsEmpty (c)))
        {
            active_set[i] = NULL;
        }
        else
        {
            failIfTrue (!mxIsDouble(c), "Matrix 'c' must be of 'double' type.");
            checkInputMatrixSize(c, num_constr[i], 1, "c");
            active_set[i] = c;
        }

    // form Ab = [A, b]
        mxArray *cat_input[4];
        mxArray *cat_output[1];

        mxArray *dimension = mxCreateDoubleMatrix (1, 1, mxREAL);
        *mxGetPr(dimension) = 2;

        cat_input[0] = dimension;
        cat_input[1] = A;
        cat_input[2] = lb;
        cat_input[3] = ub;

        if(mexCallMATLAB (1, cat_output, 4, cat_input, "cat") != 0)
        {
            mexErrMsgTxt("catenation of A and b failed!");
        }
        mxDestroyArray(dimension);

    // remember objective
        constraints[i] = cat_output[0];

        obj_type[i] = LexLS::DEFAULT_OBJECTIVE;
    }



    LexLS::TerminationStatus status;
    LexLS::LexLSI lexlsi(num_var, num_obj, num_constr.data(), obj_type.data());


// initialize solver
    try
    {
        // tolerance
        if (tolerances_are_given)
        {
            lexlsi.setTolerance(    linear_dependence_tolerance,
                                    deactivation_tolerance);
        }

        // set maximum number of iterations
        if (max_iterations_set)
        {
            lexlsi.setMaxIter(max_iterations);
        }

        // cycling handling
        if (cycling_option_enabled)
        {
            lexlsi.setCyclingHandling(max_stall_iterations);
        }

        // constraints
        if (simple_bounds_enabled)
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
            mxArray * wi;

            try
            {
                LexLS::dVectorType& w = lexlsi.getResidual(i);
                wi = mxCreateDoubleMatrix(num_constr[i], 1, mxREAL);
                for (int j = 0; j < num_constr[i]; ++j)
                {
                    mxGetPr(wi)[j] = w(j);
                }
            }
            catch (std::exception &e)
            {
                mexErrMsgTxt(e.what());
            }

            mxSetCell(output[2], i, wi);
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
            for (int j = 0; j < lexlsi_active_constraints.size(); ++j)
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


