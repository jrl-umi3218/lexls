// Time-stamp: <2014-12-09 10:12:10 drdv>
#ifndef TYPEDEFS
#define TYPEDEFS

#include <iostream>
#include <vector>
#include <algorithm> // for std::copy
#include <iterator>  // for std::ostream_iterator
#include <fstream>
#include <Eigen/Dense>

namespace LexLS
{

    /// @todo Why not unsigned?
    typedef int Index;
    typedef double RealScalar;

    /// @todo Why MatrixType, but dVectortype?
    typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixType;    
    //typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType; 
    
    typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, 1>              dVectorType;
    typedef Eigen::Matrix<RealScalar, 1, Eigen::Dynamic>              dRowVectorType;

    typedef Eigen::Matrix<     Index, Eigen::Dynamic, 1>              iVectorType;
    
    typedef Eigen::Block<MatrixType , Eigen::Dynamic, Eigen::Dynamic> dBlockType;
    typedef Eigen::Block<MatrixType , Eigen::Dynamic, 1>              dBlockType2Vector;
    typedef Eigen::VectorBlock<dVectorType, Eigen::Dynamic>           dVectorBlockType;


    enum RegularizationType
    {
        REGULARIZATION_NONE = 0,    // 0
        REGULARIZATION_TIKHONOV,    // 1
        REGULARIZATION_TIKHONOV_CG, // 2
        REGULARIZATION_R,           // 3
        REGULARIZATION_R_NO_Z,      // 4
        REGULARIZATION_RT_NO_Z,     // 5
        REGULARIZATION_RT_NO_Z_CG,  // 6
        REGULARIZATION_TIKHONOV_1,  // 7
        REGULARIZATION_TIKHONOV_2,  // 8
        REGULARIZATION_TEST         // 9
    };
    
    /**
       \brief Termination status
    */
    enum TerminationStatus
    {
        TERMINATION_STATUS_UNKNOWN = -1,
        PROBLEM_SOLVED,                  // 0
        PROBLEM_SOLVED_CYCLING_HANDLING, // 1
        MAX_NUMBER_OF_ITERATIONS_EXCEDED // 2
    };

    /**
       \brief Type of objective function 
       
       \verbatim
       Some other options (to implement ...):
       -----------------------------------------------------
       EQUALITY_OBJECTIVE           //       A*x - w  = b
       ONLY_UPPER_BOUNDS_OBJECTIVE  //       A*x - w <= bu
       ONLY_LOWER_BOUNDS_OBJECTIVE  // bl <= A*x - w
       \endverbatim

       \note Problems in ASCII files assume that DEFAULT_OBJECTIVE = 0 and SIMPLE_BOUNDS_OBJECTIVE = 1
    */
    enum ObjectiveType
    {
        DEFAULT_OBJECTIVE,         // bl <= A*x - w <= bu, assuming bl <= bu
        SIMPLE_BOUNDS_OBJECTIVE,   // bl <=   x - w <= bu, assuming bl <= bu
        SIMPLE_BOUNDS_OBJECTIVE_HP // bl <=   x     <= bu, assuming bl <= bu (only for the objective with highest priority)
    };

    /**
       \brief Type of a constraint in the working set (to be set when checking for blocking
       constraints). Used when computing the most negative/positive Lagrange multiplier.
    */
    enum ConstraintType
    {
        CONSTRAINT_TYPE_UNKNOWN, // used for initialization purposes
        LOWER_BOUND,             // bl <= A*x - w
        UPPER_BOUND,             //       A*x - w <= bu
        EQUALITY_CONSTRAINT,     //       A*x - w  = b
        CORRECT_SIGN_OF_LAMBDA   // positive if UPPER_BOUND
                                 // negative if LOWER_BOUND
    };

    /**
       \brief Performed operation during the current step
    */
    enum OperationType
    {
        UNDEFINED,         // used for initialization purposes
        ADD_CONSTRAINT,    // when constraint is added
        REMOVE_CONSTRAINT  // when constraint is removed
    };
}

#endif // TYPEDEFS
