/*
 * Copyright 2013-2021 INRIA
 */

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
    typedef unsigned int Index;
    typedef double RealScalar;

    typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> dMatrixType; // Eigen::RowMajor

    typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> dVectorType;
    typedef Eigen::Matrix<RealScalar, 1, Eigen::Dynamic> dRowVectorType;
    typedef Eigen::Matrix<     Index, Eigen::Dynamic, 1> iVectorType;

    typedef Eigen::Block<dMatrixType , Eigen::Dynamic, Eigen::Dynamic> dBlockType;
    typedef Eigen::Block<dMatrixType , Eigen::Dynamic, 1>              dBlockType2Vector;
    typedef Eigen::VectorBlock<dVectorType, Eigen::Dynamic>            dVectorBlockType;


    enum RegularizationType
    {
        REGULARIZATION_NONE = 0,       // 0
        REGULARIZATION_TIKHONOV,       // 1
        REGULARIZATION_TIKHONOV_CG,    // 2
        REGULARIZATION_R,              // 3
        REGULARIZATION_R_NO_Z,         // 4
        REGULARIZATION_RT_NO_Z,        // 5
        REGULARIZATION_RT_NO_Z_CG,     // 6
        REGULARIZATION_TIKHONOV_1,     // 7
        REGULARIZATION_TIKHONOV_2,     // 8
        REGULARIZATION_TEST            // 9
    };

    /**
       \brief Termination status
    */
    enum TerminationStatus
    {
        TERMINATION_STATUS_UNKNOWN = -1,     // -1
        PROBLEM_SOLVED,                      // 0
        PROBLEM_SOLVED_CYCLING_HANDLING,     // 1
        MAX_NUMBER_OF_FACTORIZATIONS_EXCEEDED // 2
    };

    /**
       \brief Type of objective function
    */
    enum ObjectiveType
    {
        GENERAL_OBJECTIVE = 0,    // bl <= A*x - w <= bu, assuming bl <= bu
        SIMPLE_BOUNDS_OBJECTIVE   // bl <=   x - w <= bu, assuming bl <= bu (only for the objective with highest priority)
    };

    /**
       \brief Type of a constraint in the working set.
    */
    enum ConstraintActivationType
    {
        CTR_INACTIVE = 0,        // 0
        CTR_ACTIVE_LB,           // 1: bl <= A*x - w
        CTR_ACTIVE_UB,           // 2:       A*x - w <= bu
        CTR_ACTIVE_EQ,           // 3:       A*x - w  = b
        CORRECT_SIGN_OF_LAMBDA   // 4: positive if CTR_ACTIVE_UB, negative if CTR_ACTIVE_LB (for internal use)
    };

    class ParametersLexLSE
    {
    public:
        /**
            \brief Tolerance: linear dependence (used when solving an LexLSE problem)
        */
        RealScalar tol_linear_dependence;

        /**
            \brief Max number of iterations for cg_tikhonov(...)

            \note used only with regularization_type = REGULARIZATION_TIKHONOV_CG
        */
        Index max_number_of_CG_iterations;

        /**
            \brief Type of regularization (Tikhonov, Basic Tikhonov, ...)
        */
        RegularizationType regularization_type;

        /**
         * @brief
         * @todo add documentation
         */
        RealScalar variable_regularization_factor;

        ParametersLexLSE()
        {
            setDefaults();
        }

        void print()
        {
            printf("tol_linear_dependence          = %e \n", tol_linear_dependence);
            printf("regularization_type            = %d \n", regularization_type);
            printf("variable_regularization_factor = %e \n", variable_regularization_factor);
            printf("max_number_of_CG_iterations    = %d \n", max_number_of_CG_iterations);
            printf("\n");
        }

        void setDefaults()
        {
            tol_linear_dependence          = 1e-12;
            max_number_of_CG_iterations    = 10;
            regularization_type            = REGULARIZATION_NONE;
            variable_regularization_factor = 0.0;
        }
    };

    class ParametersLexLSI
    {
    public:
        /**
           \brief Maximum number of factorizations (if reached, the solver terminates)
        */
        Index max_number_of_factorizations;

        /**
            \brief Tolerance: linear dependence (used when solving an LexLSE problem)
        */
        RealScalar tol_linear_dependence;

        /**
            \brief Tolerance: absolute value of Lagrange multiplier to be considered with "wrong" sign
        */
        RealScalar tol_wrong_sign_lambda;

        /**
            \brief Tolerance: absolute value of Lagrange multiplier to be considered with "correct" sign
        */
        RealScalar tol_correct_sign_lambda;

        /**
            \brief Tolerance: used to determine whether a constraint has been violated

            \note This tolerance is used when checking for blocking constraint and when initializing v0.
        */
        RealScalar tol_feasibility;

        /**
            \brief Type of regularization (Tikhonov, Basic Tikhonov, ...)
        */
        RegularizationType regularization_type;

        /**
            \brief Max number of iterations for cg_tikhonov(...)

            \note used only with regularization_type = REGULARIZATION_TIKHONOV_CG
        */
        Index max_number_of_CG_iterations;

        /**
         * \brief When variable_regularization_factor = 0 the user specified regularization factors
         * are used directly. When variable_regularization_factor != 0 an estimation of the
         * conditioning (conditioning_estimate) of each level is made (during the LOD). Then if
         * conditioning_estimate > variable_regularization_factor no regularization is applied,
         * while if conditioning_estimate < variable_regularization_factor the regularization
         * factors (provided by the user) are modified (there are various approaches to do this, see
         * lexlse::factorize(...)).
         *
         * \attention This functionality is not mature yet (use with caution).
         */
        RealScalar variable_regularization_factor;

        /**
            \brief If cycling_handling_enabled == true, cycling handling is performed
        */
        bool cycling_handling_enabled;

        /**
         * \brief Maximum number of attempts for cycling handling
         *
         * @todo cycling is not always detected
         */
        Index cycling_max_counter;

        /**
         * \brief Ammount of relaxation performed during each attempt to handle cycling
         */
        RealScalar cycling_relax_step;

        /**
         * \brief Used to output information for intermediate iterations of the solver
         */
        std::string output_file_name;

        /**
            \brief Allows modification of the user guess for x (see doc/hot_start.pdf)
        */
        bool modify_x_guess_enabled;

        /**
            \brief Allows modification of the user guess for active constraints (see doc/hot_start.pdf)
        */
        bool modify_type_active_enabled;

        /**
            \brief Allows modification of the user guess for inactive constraints (see doc/hot_start.pdf)
        */
        bool modify_type_inactive_enabled;

        /**
            \brief Generate the smallest possible v0 (see doc/hot_start.pdf)
        */
        bool set_min_init_ctr_violation;

        /**
            \brief If true, use phase1_v0() instead of phase1()
        */
        bool use_phase1_v0;

        /**
            \brief If true, gather information about activations and deactivations
        */
        bool log_working_set_enabled;

        /**
           \brief If true, deactivate first constraints with lambda with wrong sign. Otherwise, deactivate
           constraints with largest lambda (with wrong sign).
        */
        bool deactivate_first_wrong_sign;

        ParametersLexLSI()
        {
            setDefaults();
        }

        void print()
        {
            printf("max_number_of_factorizations   = %d \n", max_number_of_factorizations);
            printf("tol_linear_dependence          = %e \n", tol_linear_dependence);
            printf("tol_wrong_sign_lambda          = %e \n", tol_wrong_sign_lambda);
            printf("tol_correct_sign_lambda        = %e \n", tol_correct_sign_lambda);
            printf("tol_feasibility                = %e \n", tol_feasibility);
            printf("cycling_handling_enabled       = %d \n", cycling_handling_enabled);
            printf("cycling_max_counter            = %d \n", cycling_max_counter);
            printf("cycling_relax_step             = %e \n", cycling_relax_step);
            printf("regularization_type            = %d \n", regularization_type);
            printf("max_number_of_CG_iterations    = %d \n", max_number_of_CG_iterations);
            printf("variable_regularization_factor = %e \n", variable_regularization_factor);
            printf("modify_x_guess_enabled         = %d \n", modify_x_guess_enabled);
            printf("modify_type_active_enabled     = %d \n", modify_type_active_enabled);
            printf("modify_type_inactive_enabled   = %d \n", modify_type_inactive_enabled);
            printf("set_min_init_ctr_violation     = %d \n", set_min_init_ctr_violation);
            printf("use_phase1_v0                  = %d \n", use_phase1_v0);
            printf("log_working_set_enabled        = %d \n", log_working_set_enabled);
            printf("deactivate_first_wrong_sign    = %d \n", deactivate_first_wrong_sign);
            printf("\n");
        }

        void setDefaults()
        {
            max_number_of_factorizations   = 200;

            tol_linear_dependence          = 1e-12;
            tol_wrong_sign_lambda          = 1e-08;
            tol_correct_sign_lambda        = 1e-12;
            tol_feasibility                = 1e-13;

            cycling_handling_enabled       = false;
            cycling_max_counter            = 50;
            cycling_relax_step             = 1e-08;

            regularization_type            = REGULARIZATION_NONE;
            max_number_of_CG_iterations    = 10;
            variable_regularization_factor = 0.0;

            modify_x_guess_enabled         = false;
            modify_type_active_enabled     = false;
            modify_type_inactive_enabled   = false;
            set_min_init_ctr_violation     = true;

            use_phase1_v0                  = false;
            log_working_set_enabled        = false;

            deactivate_first_wrong_sign    = false;
        }
    };

    /**
       \brief A class for handling exceptions
    */
    class Exception: public std::exception
    {
    public:

        explicit Exception(const char *message): ExceptionMessage(message) {}

        ~Exception() throw() {}

        const char* what() const throw()
        {
            return ExceptionMessage.c_str();
        }

    private:

        std::string ExceptionMessage;
    };

    /**
       \brief A class used to store information related to constraints.

       \todo This class is very similar to WorkingSetLogEntry and ConstraintIdentifier (to organize
       better). Now, I have this separate class just for testing purposes.
    */

    class ConstraintInfo
    {
    public:
        ConstraintInfo():
            obj_index(0),
            ctr_index(0){}

        ConstraintInfo(int obj_index_, int ctr_index_)
        {
            obj_index = obj_index_;
            ctr_index = ctr_index_;
        }

        void increment_obj_index(int increment)
        {
            obj_index += increment;
        }

        void set_ctr_index(int ctr_index_)
        {
            ctr_index = ctr_index_;
        }

        int get_obj_index()
        {
            return obj_index;
        }

        int get_ctr_index()
        {
            return ctr_index;
        }

        void print() const
        {
            printf("obj_index = %d, ctr_index = %d \n", obj_index, ctr_index);
        }

        friend bool operator== (const ConstraintInfo &c1, const ConstraintInfo &c2);

    private:
        int obj_index;
        int ctr_index;
    };

    bool operator== (const ConstraintInfo &c1, const ConstraintInfo &c2)
    {
        bool isEqual = false;
        if ((c1.obj_index == c2.obj_index) &&
            (c1.ctr_index == c2.ctr_index))
        {
            isEqual = true;
        }

        return isEqual;
    }


    /**
        \brief A class used to store information related to the working set.
    */
    class WorkingSetLogEntry
    {
    public:

        WorkingSetLogEntry(){}

        WorkingSetLogEntry(Index obj_index_, Index ctr_index_, ConstraintActivationType ctr_type_,
                           RealScalar alpha_or_lambda_, Index rank_):
            obj_index(obj_index_),
            ctr_index(ctr_index_),
            ctr_type(ctr_type_),
            alpha_or_lambda(alpha_or_lambda_),
            rank(rank_),
            cycling_detected(false){}

        /**
            \brief Index of objective
        */
        Index obj_index;

        /**
            \brief Index of constraint

            \note This index could mean different things (depending on how an instance of this class
            is used): (i) if a constraint is to be included in the active set ctr_index is the index
            within objective obj_index; (ii) if a constraint is to be removed from the working set,
            ctr_index indicates the index within the set of active constraints.
        */
        Index ctr_index;

        /**
            \brief Type of constraint
        */
        ConstraintActivationType ctr_type;

        /**
           \brief Used to store the step length when constraint is added and largest (in
           absolute value) worng lambda when constraint is removed

           \note: not used when comparing
        */
        RealScalar alpha_or_lambda;

        /**
           \brief Rank of active constraints
        */
        Index rank;

        /**
           \brief it true, cycling has been detected
        */
        bool cycling_detected;
    };


    /**
        \brief A class used to identify a constraint.

        \todo Remove stuff not used in this class (alpha_or_lambda, cycling_detected)
    */
    class ConstraintIdentifier
    {
    public:

        ConstraintIdentifier(){}

        ConstraintIdentifier(Index obj_index_, Index ctr_index_, ConstraintActivationType ctr_type_):
            obj_index(obj_index_),
            ctr_index(ctr_index_),
            ctr_type(ctr_type_) {}

        ConstraintIdentifier(Index obj_index_, Index ctr_index_, ConstraintActivationType ctr_type_,
                             RealScalar alpha_or_lambda_):
            obj_index(obj_index_),
            ctr_index(ctr_index_),
            ctr_type(ctr_type_),
            alpha_or_lambda(alpha_or_lambda_),
            cycling_detected(false) {}

        void set(Index obj_index_, Index ctr_index_, ConstraintActivationType ctr_type_)
        {
            obj_index = obj_index_;
            ctr_index = ctr_index_;
            ctr_type  = ctr_type_;
        }

        bool operator==(const ConstraintIdentifier& ci)
        {
            return compare(ci);
        }

        bool compare(const ConstraintIdentifier& ci)
        {
            if (obj_index != ci.obj_index)
            {
                return false;
            }

            if (ctr_index != ci.ctr_index)
            {
                return false;
            }

            if (ctr_type != ci.ctr_type)
            {
                return false;
            }

            return true;
        }

        void print()
        {
            printf("type = %d, obj = %2d, ind = %3d \n", ctr_type, obj_index, ctr_index);
        }


        /**
            \brief Index of objective
        */
        Index obj_index;

        /**
            \brief Index of constraint

            \note This index could mean different things (depending on how an instance of this class
            is used): (i) if a constraint is to be included in the active set ctr_index is the index
            within objective obj_index; (ii) if a constraint is to be removed from the working set,
            ctr_index indicates the index within the set of active constraints.
        */
        Index ctr_index;

        /**
            \brief Type of constraint
        */
        ConstraintActivationType ctr_type;

        /**
           \brief Used to store the step length when constraint is added and largest (in
           absolute value) worng lambda when constraint is removed

           \note: not used when comparing
        */
        RealScalar alpha_or_lambda;

        /**
           \brief it true, cycling has been detected
        */
        bool cycling_detected;
    };


    // ----------------------------------------------------------------------------------------------------------
    // internal
    // ----------------------------------------------------------------------------------------------------------
    namespace internal
    {
        /**
           \brief Performed operation during the current step
        */
        enum OperationType
        {
            OPERATION_UNDEFINED, // used for initialization purposes
            OPERATION_ADD,       // when constraint is added
            OPERATION_REMOVE     // when constraint is removed
        };

        /**
            \brief A single Given's rotation (I use this class as a structure)
        */
        class GivensRotation
        {
        public:
            GivensRotation(){}

            GivensRotation(RealScalar a, RealScalar b, Index i_, Index j_)
            {
                set(a, b, i_, j_);
            }

            void set(RealScalar a, RealScalar b, Index i_, Index j_)
            {
                G.makeGivens(a,b);
                i = i_;
                j = j_;
            }

            void print()
            {
                printf("i = % d, j = % d, c = % f, s = % f \n", i, j, G.c(), G.s());
            }

            Eigen::JacobiRotation<RealScalar> G;
            Index i;
            Index j;
        };

        /**
            \brief A sequence of Given's rotations
        */
        class GivensRotationSequence
        {
        public:
            GivensRotationSequence(){}

            GivensRotationSequence(Index dim)
            {
                seq.reserve(dim);
            }

            void push(GivensRotation& G)
            {
                seq.push_back(G);
            }

            Eigen::JacobiRotation<RealScalar>& get(Index k)
            {
                return seq[k].G;
            }

            Index get_i(Index k)
            {
                return seq[k].i;
            }

            Index get_j(Index k)
            {
                return seq[k].j;
            }

            Index size()
            {
                return seq.size();
            }

        private:
            std::vector<GivensRotation> seq;
        };

        /**
            \brief Information about an objective in a LexLSE problem

            \note I use this as a data-structure and all fields are public.
        */
        class ObjectiveInfo
        {
        public:

            ObjectiveInfo():
                dim(0),
                rank(0),
                first_row_index(0),
                first_col_index(0),
                regularization_factor(0.0){}

            /**
               \brief Print objective information.
            */
            void print() const
            {
                printf("first_row_index = %d, first_col_index = %d, dim = %d, rank = %d, regularization_factor = %f \n", first_row_index, first_col_index, dim, rank, regularization_factor);
            }

            /*
              \brief Number of constraints involved in (LexLSE) objective

              \note set during initialization.
            */
            Index dim;

            /*
              \brief Rank of constraints involved in (LexLSE) objective

              \attention The rank is estimated during the factorization step in a naive way (and could
              be wrong) - see p. 260, Example 5.5.1 of "Matrix computations" by Golub & van Loan.
            */
            Index rank;

            /*
              \brief Initial row index - depends only on the number of constraints involved in (LexLSE) objectives

              \note computed during initialization.
            */
            Index first_row_index;

            /*
              \brief Initial column index - depends on the ranks of constraints involved in (LexLSE) objectives

              \note computed during factorization.
            */
            Index first_col_index;

            /*
              \brief Regularization factor for the current objective (default: 0.0)
            */
            RealScalar regularization_factor;
        };

    } // END namespace internal
    // ----------------------------------------------------------------------------------------------------------

} // END namespace LexLS

#endif // TYPEDEFS
