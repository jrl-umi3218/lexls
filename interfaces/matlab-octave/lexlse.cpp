/*
 * Copyright 2013-2021 INRIA
 */

#include <cmath>
#include <cstdarg>
#include <vector>

#include <lexls/lexlse.h>
#include <lexls/typedefs.h>
#include <mex.h>
#include <stdint.h>

#include "lexls_common.h"

const int MIN_INPUTS                  = 1;
const int MAX_INPUTS                  = 2;
const int MIN_OUTPUTS                 = 1;
const int MAX_OUTPUTS                 = 3;
const int MIN_NUMBER_OF_FIELDS_IN_OBJ = 2;

/**
 * @brief The main function.
 */
void mexFunction(int num_output, mxArray *output[], int num_input, const mxArray *input[])
{
    checkInputOutput(num_output, output, num_input, input, MIN_OUTPUTS, MAX_OUTPUTS, MIN_INPUTS, MAX_INPUTS,
                     MIN_NUMBER_OF_FIELDS_IN_OBJ);

    // parameters of the solver
    LexLS::ParametersLexLSE lexlse_parameters;
    lexlse_parameters.setDefaults();

    bool is_variables_fixing_enabled = false;

    unsigned int get_least_norm_solution = 0;

    bool is_regularization_set = false;
    std::vector<double> regularization_factors;

    // parse parameters

    if (num_input == 2)
    {
        const mxArray *options_struct = input[1];

        // ================================================
        // parse options

        getOptionDouble(&lexlse_parameters.tol_linear_dependence, options_struct, "tol_linear_dependence");

        getOptionBool(&is_variables_fixing_enabled, options_struct, "enable_fixed_variables");

        is_regularization_set = getOptionArray(regularization_factors, mxGetNumberOfElements(input[0]), options_struct,
                                               "regularization_factors");

        getOptionUnsignedInteger(&get_least_norm_solution, options_struct, "get_least_norm_solution");

        unsigned int regularization_type = 0;
        if (getOptionUnsignedInteger(&regularization_type, options_struct, "regularization_type"))
        {
            lexlse_parameters.regularization_type = static_cast<LexLS::RegularizationType>(regularization_type);
        }

        getOptionUnsignedInteger(&lexlse_parameters.max_number_of_CG_iterations, options_struct,
                                 "max_number_of_CG_iterations");

        getOptionDouble(&lexlse_parameters.variable_regularization_factor, options_struct,
                        "variable_regularization_factor");

        // ================================================
        // check provided options

        /// @todo This check should probably go to the solver interface
        if (((lexlse_parameters.regularization_type == LexLS::REGULARIZATION_NONE) && (is_regularization_set))
            || ((lexlse_parameters.regularization_type != LexLS::REGULARIZATION_NONE) && (!is_regularization_set)))
        {
            mexWarnMsgTxt("Both regularization type and regularization factors must be specified.");
        }
    }

    // parse objectives

    const mxArray *objectives = input[0];

    LexLS::Index num_obj                     = mxGetNumberOfElements(objectives);
    LexLS::Index number_of_normal_objectives = num_obj;
    LexLS::Index index_first_normal_obj      = 0;

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

        index_first_normal_obj = 1;
        --number_of_normal_objectives;

        failIfTrue(num_obj == 1, "Problems consisting of one level of fixed variables are not supported.");
    }

    LexLS::Index num_var = 0;
    int total_num_constr = 0;

    std::vector<mxArray *> constraints;
    std::vector<LexLS::Index> num_constr;
    num_constr.resize(number_of_normal_objectives);
    constraints.resize(number_of_normal_objectives);

    for (LexLS::Index i = index_first_normal_obj, index_obj = 0; index_obj < number_of_normal_objectives;
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
        // remember objective
        constraints[index_obj] = catenateMatrices(A, b);
    }

    // initialize solver
    try
    {
        LexLS::internal::LexLSE lexlse(num_var, number_of_normal_objectives, num_constr.data());
        lexlse.setParameters(lexlse_parameters);

        // fixed variables
        if (is_variables_fixing_enabled)
        {
            int fixed_var_num = mxGetM(fixed_var_i);

            lexlse.setFixedVariablesCount(fixed_var_num);
            for (int i = 0; i < fixed_var_num; ++i)
            {
                LexLS::Index index = static_cast<int>(round(mxGetPr(fixed_var_i)[i])) - 1;
                //if ((index < 0) || (index > num_var))
                if (index > num_var)
                {
                    mexErrMsgTxt("Index of a fixed variable is out of bounds!");
                }
                lexlse.fixVariable(index, mxGetPr(fixed_var_b)[i]);
            }
        }

        if (is_regularization_set)
        {
            for (LexLS::Index i = index_first_normal_obj, j = 0; i < num_obj; ++i, ++j)
            {
                lexlse.setRegularizationFactor(j, regularization_factors[i]);
            }
        }

        // constraints
        for (LexLS::Index i = 0; i < number_of_normal_objectives; ++i)
        {
            lexlse.setData(i, Eigen::Map<LexLS::dMatrixType>(mxGetPr(constraints[i]), num_constr[i], num_var + 1));
        }

        // solve the problem
        lexlse.factorize();

        if (get_least_norm_solution == 1)
        {
            lexlse.solveLeastNorm_1();
        }
        else if (get_least_norm_solution == 2)
        {
            lexlse.solveLeastNorm_2();
        }
        else if (get_least_norm_solution == 3)
        {
            lexlse.solveLeastNorm_3();
        }
        else
        {
            lexlse.solve();
        }

        // output solution
        LexLS::dVectorType &x = lexlse.get_x();
        output[0]             = mxCreateDoubleMatrix(num_var, 1, mxREAL);
        double *x_out         = mxGetPr(output[0]);
        for (LexLS::Index i = 0; i < num_var; ++i)
        {
            x_out[i] = x(i);
        }

        // output additional information
        if (num_output >= 2)
        {
            int num_info_fields            = 1;
            const char *info_field_names[] = {"status"};

            output[1] = mxCreateStructMatrix(1, 1, num_info_fields, info_field_names);

            mxArray *info_status                   = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            ((INT32_T *)mxGetData(info_status))[0] = STATUS_OK;

            mxSetField(output[1], 0, "status", info_status);
        }

        // output residual
        if (num_output >= 3)
        {
            LexLS::dVectorType &w = lexlse.get_v();

            output[2]   = mxCreateCellMatrix(number_of_normal_objectives, 1);
            int index_w = 0;
            for (LexLS::Index i = 0; i < number_of_normal_objectives; ++i)
            {
                mxArray *wi = mxCreateDoubleMatrix(num_constr[i], 1, mxREAL);
                for (LexLS::Index j = 0; j < num_constr[i]; ++j)
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
        mexErrMsgTxt(e.what());
    }
}
