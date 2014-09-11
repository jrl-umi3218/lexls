/**
    @file
    @author  Alexander Sherikov

    @brief 
*/


const int MIN_INPUTS = 1;
const int MAX_INPUTS = 2;

enum ReturnStatus
{
    STATUS_OK = 0,
    STATUS_CYCLING_HANDLING,
    STATUS_MAXITER,
    STATUS_UNKNOWN_FAILURE
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


void checkInputMatrix(const mxArray *M, const char *name)
{
    failIfTrue (M == NULL, (std::string("Matrix '") + name + "' is missing.").c_str());
    failIfTrue (mxIsEmpty (M), (std::string("Matrix '") + name + "' is empty.").c_str());
    failIfTrue (!mxIsDouble(M), (std::string("Matrix '") + name + "' must be of 'double' type.").c_str());
}


void checkInputMatrixSize(const mxArray *M, const int nrows, const int ncols, const char *name)
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
                        const int min_number_of_fields)
{
    if ((num_input < MIN_INPUTS) || (num_input > MAX_INPUTS))
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
        if ((num_input == 2) && (!mxIsStruct(input[1])))
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
