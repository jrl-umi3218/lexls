#include <vector>
#include <cmath>
#include <cstdarg>

#include <mex.h>
#include <typedefs.h>
#include <lexlse.h>


#include "lexls_common.h"


const int MIN_INPUTS = 1;
const int MAX_INPUTS = 2;
const int MIN_OUTPUTS = 1;
const int MAX_OUTPUTS = 3;
const int MIN_NUMBER_OF_FIELDS_IN_OBJ = 2;


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
    double linear_dependence_tolerance = 0.0;
    bool is_linear_dependence_tolerance_set = false;
    bool is_variables_fixing_enabled = false;


// parse parameters

    if (num_input == 2)
    {
        const mxArray *options = input[1];

        is_linear_dependence_tolerance_set = getOptionDouble(   &linear_dependence_tolerance, 
                                                                options, 
                                                                "linear_dependence_tolerance");
        
        is_variables_fixing_enabled = getOptionBool(options, "enable_fixed_variables", is_variables_fixing_enabled);
    }

// parse objectives

    const mxArray *objectives = input[0];

    LexLS::Index num_obj = mxGetNumberOfElements (objectives);
    int index_first_normal_obj = 0;

    mxArray *fixed_var_i = NULL;
    mxArray *fixed_var_b = NULL;


    if (is_variables_fixing_enabled)
    { // fixed variables
        mxArray *A = getObjectiveMatrix(objectives, 0, "A");
        mxArray *b = getObjectiveMatrix(objectives, 0, "b");
            

    // check A and b
        int num_rows = mxGetM(A);
        checkInputMatrixSize(A, num_rows, 1, 0, "A");
        checkInputMatrixSize(b, num_rows, 1, 0, "b");

        fixed_var_i = A;
        fixed_var_b = b;

        is_variables_fixing_enabled = true;
        index_first_normal_obj = 1;
        --num_obj;
    }


    LexLS::Index num_var = 0;
    int total_num_constr = 0;

    std::vector<mxArray *> constraints;
    std::vector<LexLS::Index> num_constr;
    num_constr.resize(num_obj);
    constraints.resize(num_obj);


    for (int i = index_first_normal_obj, index_obj = 0; 
            index_obj < num_obj; 
            ++index_obj, ++i)
    {
        mxArray *A = getObjectiveMatrix(objectives, i, "A");
        mxArray *b = getObjectiveMatrix(objectives, i, "b");

    // check A and b
        num_constr[index_obj] = mxGetM(A);
        total_num_constr += num_constr[index_obj];

        if (i == index_first_normal_obj)
        {
            num_var = mxGetN(A);
        }

        checkInputMatrixSize(A, num_constr[index_obj], num_var, i, "A");
        checkInputMatrixSize(b, num_constr[index_obj], 1, i, "b");


    // form Ab = [A, b]
        mxArray *cat_input[3];
        mxArray *cat_output[1];

        mxArray *dimension = mxCreateDoubleMatrix (1, 1, mxREAL);
        *mxGetPr(dimension) = 2;

        cat_input[0] = dimension;
        cat_input[1] = A;
        cat_input[2] = b;

        if(mexCallMATLAB (1, cat_output, 3, cat_input, "cat") != 0)
        {
            mexErrMsgTxt("Catenation of A and b failed!");
        }
        mxDestroyArray(dimension);

    // remember objective
        constraints[index_obj] = cat_output[0];
    }


// initialize solver
    try
    { 
        LexLS::LexLSE lexlse(num_var, num_obj, num_constr.data());

        // tolerance
        if (is_linear_dependence_tolerance_set)
        {
            lexlse.setTolerance(linear_dependence_tolerance);
        }

        // fixed variables
        if (is_variables_fixing_enabled)
        {
            int fixed_var_num = mxGetM(fixed_var_i);

            lexlse.setFixedVariablesCount(fixed_var_num);
            for (int i = 0; i < fixed_var_num; ++i)
            {
                LexLS::Index index = static_cast <int> (round(mxGetPr(fixed_var_i)[i])) - 1;
                if ((index < 0) || (index > num_var))
                {
                    mexErrMsgTxt("Index of a fixed variable is out of bounds!");
                }
                lexlse.fixVariable(index, mxGetPr(fixed_var_b)[i]);
            }
        }


        // constraints
        for (int i = 0; i < num_obj; ++i)
        {
            lexlse.setData(
                    i, 
                    Eigen::Map<LexLS::MatrixType>(
                        mxGetPr(constraints[i]), 
                        num_constr[i], 
                        num_var + 1));
        }


    // solve the problem
        lexlse.factorize();
        lexlse.solve();
    

    // output solution
        LexLS::dVectorType& x = lexlse.get_x();
        output[0] = mxCreateDoubleMatrix(num_var, 1, mxREAL);
        double *x_out = mxGetPr(output[0]);
        for (int i = 0; i < num_var; ++i)
        {
            x_out[i] = x(i);
        }

    // output additional information
        if (num_output >= 2)
        {
            int num_info_fields = 1; 
            const char *info_field_names[] = {"status"}; 

            output[1] = mxCreateStructMatrix(1, 1, num_info_fields, info_field_names);

            mxArray *info_status = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            ((INT32_T *) mxGetData (info_status))[0] = STATUS_OK;

            mxSetField (output[1], 0, "status", info_status);
        }

    // output residual
        if (num_output >= 3)
        {
            LexLS::dVectorType& w = lexlse.getResidual();

            output[2] = mxCreateCellMatrix(num_obj, 1);
            int index_w = 0;
            for (int i = 0; i < num_obj; ++i)
            {
                mxArray * wi = mxCreateDoubleMatrix(num_constr[i], 1, mxREAL);
                for (int j = 0; j < num_constr[i]; ++j)
                {
                    mxGetPr(wi)[j] = w(index_w);
                    ++index_w;
                }
                mxSetCell(output[2], i, wi);
            }
        }
    }
    catch (std::exception &e)
    {
        mexErrMsgTxt (e.what());
    }
}
