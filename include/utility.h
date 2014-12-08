// Time-stamp: <2014-12-08 09:37:07 drdv>
#ifndef UTILITY
#define UTILITY

#include <typedefs.h>

/**
   \brief Contains routines that can be used to solve a lexicographic least-squares problem with
   inequality constrains
*/
namespace LexLS
{
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

    // ----------------------------------------------------------------------------

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
    class ObjInfoType
    {    
    public:

        ObjInfoType():
            dim(0),
            rank(0),
            FirstRowIndex(0), 
            FirstColIndex(0),
            ObjType(DEFAULT_OBJECTIVE) {}

        /**
           \brief Print objective information.
        */                                        
        void print() const
        {
            printf("FirstRowIndex = %d, FirstColIndex = %d, dim = %d, rank = %d \n", FirstRowIndex, FirstColIndex, dim, rank);
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
        Index FirstRowIndex; 
    
        /*
          \brief Initial column index - depends on the ranks of constraints involved in (LexLSE) objectives
      
          \note computed during factorization. 
        */
        Index FirstColIndex; 

        /** 
            \brief Type of the objective

            \note Currently only DEFAULT_OBJECTIVE and SIMPLE_BOUNDS_OBJECTIVE_HP are handled
        */
        ObjectiveType ObjType;
    };

    /** 
        \brief A class used to identify a constraint (used only in the cycling detection).
    */
    class ConstraintIdentifierType
    {    
    public:
    
        ConstraintIdentifierType(){}

        ConstraintIdentifierType(Index ObjIndex_, Index CtrIndex_, ConstraintType CtrType_):
            ObjIndex(ObjIndex_),
            CtrIndex(CtrIndex_),
            CtrType(CtrType_) {}
        
        void set(Index ObjIndex_, Index CtrIndex_, ConstraintType CtrType_)
        {
            ObjIndex = ObjIndex_;
            CtrIndex = CtrIndex_;
            CtrType  = CtrType_;
        }

        bool operator==(const ConstraintIdentifierType& ci)
        {
            return compare(ci);
        }

        bool compare(const ConstraintIdentifierType& ci)
        {
            if (ObjIndex != ci.ObjIndex)
                return false;

            if (CtrIndex != ci.CtrIndex)
                return false;

            if (CtrType != ci.CtrType)
                return false;

            return true;
        }

        void print()
        {
            printf("ObjIndex = %d, CtrIndex = %d, CtrType = %d\n", ObjIndex, CtrIndex, CtrType);
        }

        Index getObjIndex()
        {
            return ObjIndex;
        }

        Index getCtrIndex()
        {
            return CtrIndex;
        }

        ConstraintType getCtrType()
        {
            return CtrType;
        }

        /** 
            \brief Index of objective
        */ 
        Index ObjIndex;
        
        /** 
            \brief Index of constraint

            \note This index could mean different things (depending on how an instance of this class
            is used): (i) if a constraint is to be included in the active set CtrIndex is the index
            within objective ObjIndex; (ii) if a constraint is to be removed from the working set,
            CtrIndex indicates the index within the set of active constraints.
        */ 
        Index CtrIndex;

        /** 
            \brief Type of constraint
        */ 
        ConstraintType CtrType;
    };
    
    /**
       \brief Delete the content of an ASCII file
    */
    inline void flushFile(const char* FileName)
    {
        std::ofstream file(FileName);
    }

    /**
       \brief Compare doubles
    */
    inline bool isEqual(RealScalar a, RealScalar b, RealScalar TOLERANCE_FOR_EQUALITY_BETWEEN_DOUBLES=1e-15)
    {
        return std::abs(a - b) < TOLERANCE_FOR_EQUALITY_BETWEEN_DOUBLES;
    }
    
    /**
       \brief generate random double 
    */
    inline double rand_double(double min=-1, double max=1)
    {
        double d = (double) rand() / (double) RAND_MAX;
        return min + d*(max - min);
    }

    /**
       \brief Matlab does not display cout ... 
    */
    inline void print_eigen_matrix(MatrixType M, const char* variable_name)
    {
        printf(" %s(%d,%d) = \n",variable_name,(Index)M.rows(),(Index)M.cols());
        for (Index d1=0; d1<M.rows(); d1++)
        {
            for (Index d2=0; d2<M.cols(); d2++)
            {
                printf(" % .7f ",M(d1,d2));
            }
            printf("\n");
        }
        printf("\n");
    }
    
} // END namespace LexLS

#endif // UTILITY
