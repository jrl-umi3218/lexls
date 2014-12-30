#include <vector>
#include <cmath>
#include <cstdarg>

#include <stdint.h>
#include <mex.h>
#include <typedefs.h>
#include <lexlse.h>


#include "lexls_common.h"


const int MIN_INPUTS = 1;
const int MAX_INPUTS = 2;
const int MIN_OUTPUTS = 1;
const int MAX_OUTPUTS = 3;
const int MIN_NUMBER_OF_FIELDS_IN_OBJ = 2;


class LexlseOptions
{
    public:
        bool is_linear_dependence_tolerance_set;
        double linear_dependence_tolerance;

        bool is_variables_fixing_enabled;

        unsigned int get_least_norm_solution;

        bool is_regularization_set;
        std::vector<double> regularization;

        bool is_regularization_type_set;
        LexLS::RegularizationType regularizationType;

        bool is_regularizationMaxIterCG_set;
        unsigned int regularizationMaxIterCG;

        LexlseOptions()
        {
            is_linear_dependence_tolerance_set = false;
            linear_dependence_tolerance = 0.0;

            is_variables_fixing_enabled = false;

            is_regularization_set = false;

            is_regularization_type_set = false;

            get_least_norm_solution = 0;

            regularizationMaxIterCG = 10;

            regularizationType = LexLS::REGULARIZATION_NONE;
        }
};



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
    LexlseOptions options;


// parse parameters

    if (num_input == 2)
    {
        const mxArray *options_struct = input[1];


        // ================================================
        // parse options

        options.is_linear_dependence_tolerance_set = getOptionDouble( &options.linear_dependence_tolerance, 
                                                                      options_struct, 
                                                                      "tolLinearDependence");
        
        getOptionBool(&options.is_variables_fixing_enabled,    
                      options_struct, 
                      "enable_fixed_variables");

        options.is_regularization_set = getOptionArray( options.regularization, 
                                                        mxGetNumberOfElements (input[0]),
                                                        options_struct, 
                                                        "regularization");

        getOptionUnsignedInteger(   &options.get_least_norm_solution, 
                            options_struct, 
                            "get_least_norm_solution");

        unsigned int regularization_type = 0;
        options.is_regularization_type_set = getOptionUnsignedInteger( &regularization_type, 
                                                               options_struct, 
                                                               "regularizationType");
        if (options.is_regularization_type_set)
        {
            options.regularizationType = static_cast <LexLS::RegularizationType> (regularization_type);
        }

        options.is_regularizationMaxIterCG_set = getOptionUnsignedInteger(   &options.regularizationMaxIterCG, 
                                                                     options_struct, 
                                                                     "regularizationMaxIterCG");

        // ================================================
        // check provided options

        /// @todo This check should probably go to the solver interface
        if ( 
                ((options.regularizationType == LexLS::REGULARIZATION_NONE) && (options.is_regularization_set)) ||
                ((options.regularizationType != LexLS::REGULARIZATION_NONE) && (!options.is_regularization_set))
           )
        {
            mexWarnMsgTxt("Both regularization type and regularization factors must be specified.");
        }
    }

// parse objectives

    const mxArray *objectives = input[0];

    LexLS::Index num_obj = mxGetNumberOfElements (objectives);
    LexLS::Index number_of_normal_objectives = num_obj;
    int index_first_normal_obj = 0;

    mxArray *fixed_var_i = NULL;
    mxArray *fixed_var_b = NULL;


    if (options.is_variables_fixing_enabled)
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


    for (int i = index_first_normal_obj, index_obj = 0; 
         index_obj < number_of_normal_objectives; 
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
        LexLS::LexLSE lexlse(num_var, number_of_normal_objectives, num_constr.data());

        // tolerance
        if (options.is_linear_dependence_tolerance_set)
        {
            lexlse.setTolerance(options.linear_dependence_tolerance);
        }

        // fixed variables
        if (options.is_variables_fixing_enabled)
        {
            int fixed_var_num = mxGetM(fixed_var_i);

            lexlse.setFixedVariablesCount(fixed_var_num);
            for (int i = 0; i < fixed_var_num; ++i)
            {
                LexLS::Index index = static_cast <int> (round(mxGetPr(fixed_var_i)[i])) - 1;
                //if ((index < 0) || (index > num_var))
                if (index > num_var)
                {
                    mexErrMsgTxt("Index of a fixed variable is out of bounds!");
                }
                lexlse.fixVariable(index, mxGetPr(fixed_var_b)[i]);
            }
        }

        if (options.is_regularization_set)
        {
            for (int i = index_first_normal_obj, j = 0; i < num_obj; ++i, ++j)
            {
                lexlse.setRegularization(j, options.regularization[i]);
            }
        }

        if (options.is_regularization_type_set)
        {
            lexlse.setRegularizationType(options.regularizationType);
        }

        if (options.is_regularizationMaxIterCG_set)
        {
            lexlse.setRegularizationMaxIterCG(options.regularizationMaxIterCG);
        }

        // constraints
        for (int i = 0; i < number_of_normal_objectives; ++i)
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

        if (options.get_least_norm_solution == 1) 
        {
            lexlse.solveLeastNorm_1();
        }
        else if (options.get_least_norm_solution == 2) 
        {
            lexlse.solveLeastNorm_2();
        }        
        else if (options.get_least_norm_solution == 3) 
        {
            lexlse.solveLeastNorm_3();
        }        
        else
        {
            lexlse.solve();
        }
        
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

            output[2] = mxCreateCellMatrix(number_of_normal_objectives, 1);
            int index_w = 0;
            for (int i = 0; i < number_of_normal_objectives; ++i)
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
