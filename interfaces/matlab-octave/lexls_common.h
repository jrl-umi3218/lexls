/**
    @file
    @author  Alexander Sherikov

    @brief 
*/

#include <sstream>


enum ReturnStatus
{
    STATUS_OK = 0,
    STATUS_CYCLING_HANDLING = 1,
    STATUS_MAXITER = 2,
    STATUS_UNKNOWN_FAILURE = 3
};


/*
void printMatrix(const mxArray *M)
{
    double *Mdata = mxGetPr(M);

    int num_rows = mxGetM(M);
    int num_cols = mxGetN(M);

    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            mexPrintf(" %5f ", Mdata[j*num_rows + i]);
        }
        mexPrintf("\n");
    }
}
*/


/**
 * @brief Prints the message and terminates execution, if the given condition is true.
 */
void failIfTrue(const bool condition, const char * message)
{
    if (condition)
    {
        mexErrMsgTxt(message);
    }
}



// ======================================================================================
// Input checks
// ======================================================================================


void checkInputMatrix(const mxArray *M, const int objective_index, const char *name)
{
    std::stringstream str_objective_index;
    str_objective_index << objective_index + 1;

    failIfTrue (M == NULL, (std::string("Matrix '") + name + "' is missing on level " + str_objective_index.str() + ".").c_str());
    failIfTrue (mxIsEmpty (M), (std::string("Matrix '") + name + "' is empty on level " + str_objective_index.str() + ".").c_str());
    failIfTrue (!mxIsDouble(M), (std::string("Matrix '") + name + "' must be of 'double' type on level " + str_objective_index.str() + ".").c_str());
}


void checkInputMatrixSize(const mxArray *M, const unsigned int nrows, const unsigned int ncols, const int objective_index, const char *name)
{
    std::stringstream str_objective_index;
    str_objective_index << objective_index + 1;

    failIfTrue (mxGetM(M) != nrows, (std::string("Wrong number of rows in matrix '") + name + "' on level " + str_objective_index.str() + ".").c_str());
    failIfTrue (mxGetN(M) != ncols, (std::string("Wrong number of columns in matrix '") + name + "' on level " + str_objective_index.str() + ".").c_str());
}


void checkInputMatrixSize(const mxArray *M, const unsigned int nrows, const unsigned int ncols, const char *name)
{
    failIfTrue (mxGetM(M) != nrows, (std::string("Wrong number of rows in matrix '") + name + "'.").c_str());
    failIfTrue (mxGetN(M) != ncols, (std::string("Wrong number of columns in matrix '") + name + "'.").c_str());
}


/**
 * @brief Checks correctness of the input / output parameters, prints a message 
 *  and terminates executione whenever an error is encountered
 */
void checkInputOutput(  const int num_output, mxArray *output[],
                        const int num_input, const mxArray *input[],
                        const int min_outputs, const int max_outputs,
                        const int min_inputs, const int max_inputs,
                        const int min_number_of_fields)
{
    if ((num_input < min_inputs) || (num_input > max_inputs))
    {
        mexErrMsgTxt("Wrong number of inputs!");
    }

    if ((num_output < min_outputs) || (num_output > max_outputs))
    {
        mexErrMsgTxt("Wrong number of outputs!");
    }

    if (!mxIsStruct(input[0]))
    {
        mexErrMsgTxt("First input parameter is not a structure.");
    }
    else
    {
        if ((num_input >= 2) && (!mxIsEmpty(input[1])) && (!mxIsStruct(input[1])))
        {
            mexErrMsgTxt("Second input parameter is not a structure.");
        }
    }

    if (mxGetNumberOfFields (input[0]) < min_number_of_fields) 
    {
        mexErrMsgTxt("Wrong number of fields in the first input structure!");
    }

    if (mxGetNumberOfElements (input[0]) < 1)
    {
        mexErrMsgTxt("At least one objective must be specified!");
    }
}


mxArray * getObjectiveMatrix(const mxArray * objectives, const int objective_index, const char *objective_matrix_id)
{
    mxArray *objective_matrix = mxGetField (objectives, objective_index, objective_matrix_id);

    // check A and b
    checkInputMatrix(objective_matrix, objective_index, objective_matrix_id);

    return (objective_matrix);
}



// ======================================================================================
// Parsing of options 
// ======================================================================================


bool getOptionDouble(double *option_value, const mxArray *option_struct, const char *option_id)
{
    bool is_parsing_successful = false;


    mxArray *option = mxGetField (option_struct, 0, option_id);

    if (option != NULL) 
    {
        // if there is such field
        failIfTrue (!mxIsDouble(option), (std::string("Option '") + option_id + "' must be of 'double' type.").c_str());
        *option_value = mxGetPr(option)[0];
        is_parsing_successful = true;
    }


    return (is_parsing_successful);
}



bool getOptionString(   std::string &option_value, 
                        const mxArray *option_struct, 
                        const char *option_id,
                        const mwSize strlen)
{
    bool is_parsing_successful = false;


    mxArray *option = mxGetField (option_struct, 0, option_id);

    if ( (option != NULL) && (strlen > 0) )
    {
        // if there is such field
        failIfTrue (!mxIsChar(option), (std::string("Option '") + option_id + "' must be of 'char' type.").c_str());

        char *char_string = new char[strlen];

        if (mxGetString(option, char_string, strlen) == 0)
        {
            option_value = char_string;
            is_parsing_successful = true;
        }

        delete char_string;
    }


    return (is_parsing_successful);
}


bool getOptionBool(bool *option_value, const mxArray *option_struct, const char *option_id)
{
    bool is_parsing_successful = false;


    mxArray * option = mxGetField (option_struct, 0, option_id);

    if (option != NULL) 
    {
        failIfTrue (!mxIsDouble(option), (std::string("Flag option '") + option_id + "' must be of 'double' type.").c_str());
        if (*mxGetPr(option) != 0)
        {
            *option_value = true;
        }
        else
        {
            *option_value = false;
        }
        is_parsing_successful = true;
    }

    return (is_parsing_successful);
}


bool getOptionUnsignedInteger(unsigned int *option_value, const mxArray *option_struct, const char *option_id)
{
    bool is_parsing_successful = false;

    double option_value_double = 0.0;


    is_parsing_successful = getOptionDouble(&option_value_double, option_struct, option_id);


    if (is_parsing_successful)
    {
        *option_value = static_cast <int> (round(option_value_double));
    }

    return (is_parsing_successful);
}


bool getOptionArray(std::vector<double> &array, 
                    const unsigned int number_of_elements, 
                    const mxArray *option_struct, 
                    const char *option_id)
{
    bool is_parsing_successful = false;


    mxArray *option = mxGetField (option_struct, 0, option_id);

    if (option != NULL) 
    {
        // if there is such field
        failIfTrue (!mxIsDouble(option), (std::string("Option '") + option_id + "' must be of 'double' type.").c_str());

        failIfTrue ((mxGetNumberOfElements(option) != number_of_elements),
                    (std::string("Wrong number of elements in option '") + option_id + "'.").c_str());

        array.clear();
        for (unsigned int i = 0; i < number_of_elements; ++i)
        {
            array.push_back( (static_cast <double *> (mxGetPr(option)) )[i]);
        }

        is_parsing_successful = true;
    }


    return (is_parsing_successful);
}



// ======================================================================================
// Operations on matrices
// ======================================================================================

mxArray * catenateMatrices( mxArray *matrix1, 
                            mxArray *matrix2, 
                            mxArray *matrix3 = NULL)
{
    int number_of_matrices = 2;

    if (matrix3 != NULL)
    {
        number_of_matrices = 3;
    }


    mxArray **cat_input = new mxArray*[number_of_matrices + 1];
    mxArray *cat_output[1];

    mxArray *dimension = mxCreateDoubleMatrix (1, 1, mxREAL);
    *mxGetPr(dimension) = 2;


    cat_input[0] = dimension;
    cat_input[1] = matrix1;
    cat_input[2] = matrix2;
    if (matrix3 != NULL)
    {
        cat_input[3] = matrix3;
    }


    if(mexCallMATLAB (1, cat_output, number_of_matrices+1, cat_input, "cat") != 0)
    {
        mexErrMsgTxt("Catenation of A and b failed!");
    }
    mxDestroyArray(dimension);
    delete cat_input;

    return (cat_output[0]);
}
