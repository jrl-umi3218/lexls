#ifndef UTILITY
#define UTILITY

#include <typedefs.h>

namespace LexLS
{    
    /**
       \brief Matlab does not display cout ... 
    */
    inline void print_eigen_matrix(dMatrixType M, const char* variable_name)
    {
        printf(" %s(%lu,%lu) = \n",variable_name,M.rows(),M.cols());
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

    // ----------------------------------------------------------------------------------------------------------
    // internal
    // ----------------------------------------------------------------------------------------------------------
    namespace internal
    {        
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

    } // END namespace internal
    // ----------------------------------------------------------------------------------------------------------
    
} // END namespace LexLS

#endif // UTILITY
