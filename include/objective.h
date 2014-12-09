// Time-stamp: <2014-12-09 15:49:20 drdv>
#ifndef OBJECTIVE
#define OBJECTIVE

#include <utility.h>
#include <workingset.h>

namespace LexLS
{
    /**
       \brief Defines an objective of a LexLSI problem

       \todo Remove throw Exception("Unknown objective type") and others like it. This should have
       to be checked at the api level.
    */
    class Objective
    {
    public:

        /**
           \brief Default constructor
        */
        Objective(): 
            nCtr(0),
            w_is_initialized(false),
            ObjType(DEFAULT_OBJECTIVE) {}

        /**
           \brief Resize the objective

           \param[in] nCtr_    Number of constraints in an objective.
           \param[in] nVar_    Number of variables (without the residuals).
           \param[in] ObjType_ Type of objective
        */
        void resize(Index nCtr_, Index nVar_, ObjectiveType ObjType_)
        {
            ObjType = ObjType_;
            nCtr    = nCtr_;
            nVar    = nVar_;

            WorkingSet.resize(nCtr);

            w.resize(nCtr);
            dw.resize(nCtr);
            
            wStar.resize(nCtr);

            Ax.resize(nCtr);
            Adx.resize(nCtr);

            if (ObjType == DEFAULT_OBJECTIVE)
            {
                data.resize(nCtr,nVar + 2); // [A,LowerBounds,UpperBounds]
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP || ObjType == SIMPLE_BOUNDS_OBJECTIVE)
            {
                data.resize(nCtr,2);         // [LowerBounds,UpperBounds]
                VarIndex.resize(nCtr);
            }
            else
            {
                throw Exception("Unknown objective type");
            }

            initialize();
        }

        /**
           \brief Given x, generate a w, such that (x,w) is a feasible initial pair for the
           constraints involved in the objective 

           \note #Ax is initialized
        */
        void phase1(dVectorType &x)
        {
            Index CtrIndex, lbIndex, ubIndex;
            ConstraintType CtrType;
            
            if (ObjType == DEFAULT_OBJECTIVE)
            {
                lbIndex = nVar;
                ubIndex = nVar+1;
                
                Ax = data.leftCols(nVar)*x; // initialize Ax
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE || ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP)
            {
                lbIndex = 0;
                ubIndex = 1;
                
                for (Index k=0; k<nCtr; k++) // initialize Ax
                    Ax(k) = x(VarIndex[k]); 
            }
            else
            {
                throw Exception("Unknown objective type");
            }

            if (!w_is_initialized)
            {
                w  = Ax - 0.5*(data.col(lbIndex)+data.col(ubIndex)); 
            
                // overwrite the residual of the active constraints
                for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                {
                    CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                    CtrType  = getActiveCtrType(CtrIndexActive);

                    if (CtrType == LOWER_BOUND)
                        w.coeffRef(CtrIndex) = Ax.coeffRef(CtrIndex) - data.coeffRef(CtrIndex,lbIndex);
                    else if (CtrType == UPPER_BOUND)
                        w.coeffRef(CtrIndex) = Ax.coeffRef(CtrIndex) - data.coeffRef(CtrIndex,ubIndex);
                }
            }

            // FIXME: Temporary hack (TO RECODE THIS FUNCTION)
            // The problem with the above implementation is that, for inactive constraints, even when 
            // lb[i] <= Ax[i] <= ub[i], a nonzero w[i] can be generated. So overwrite w[i] = 0
            /*
            for (Index CtrIndex=0; CtrIndex<nCtr; CtrIndex++)
            {
                if (!isActive(CtrIndex))
                {
                    if (Ax.coeffRef(CtrIndex) >= data.coeffRef(CtrIndex,lbIndex) && 
                        Ax.coeffRef(CtrIndex) <= data.coeffRef(CtrIndex,ubIndex))
                    {
                        //printf("w[%d]: {%f <= %f <= %f} \n",CtrIndex,data.coeffRef(CtrIndex,lbIndex), Ax.coeffRef(CtrIndex), data.coeffRef(CtrIndex,ubIndex));
                        w.coeffRef(CtrIndex) = 0.0;
                    }
                }
            }
            */

        }

        /** 
            \brief Includes in the working set the constraint with index CtrIndex (and sets its
            type)

            \param[in] CtrIndex         Index of constraint in a given LexLSI objective
            \param[in] type             Type of the constraint to be included in the working set.
        */
        void activate(Index CtrIndex, ConstraintType type)
        {
            if (CtrIndex >= nCtr)
                throw Exception("CtrIndex >= nCtr");
            
            WorkingSet.activate(CtrIndex, type);
        }

        /** 
            \brief Removes from the working set the constraint with index CtrIndexActive

            \param[in] CtrIndexActive Index of constraint in the working set in a given objective,
            i.e., Obj[ObjIndex].WorkingSet.active[CtrIndexActive] will be removed.
        */
        void deactivate(Index CtrIndexActive)
        {
            if (CtrIndexActive >= getActiveCtrCount())
            {
                throw Exception("CtrIndexActive >= number of active constraints");
            }

            WorkingSet.deactivate(CtrIndexActive);
        }

        /** 
            \brief Form an LexLSE problem (using the current working set)

            \verbatim
            -----------------------------------------
            Example (handling of simple bounds)
            -----------------------------------------
                             nVar = 6
                         VarIndex = {3,5,4}
                    active bounds = {2,0}
                  inactive bounds = {1}
                             type = {LOWER_BOUND, 
                                     UPPER_BOUND}
             -----------------------------------------
                      data(0,0)  = x(3)
                      data(1,0) <= x(5) <= data(1,1)
                                   x(4)  = data(2,1)
             -----------------------------------------
                                2 = getActiveCtrIndex(0)
                                4 = getVarIndex(2)
             -----------------------------------------
                                1 = getInactiveCtrIndex(0)
                                3 = getVarIndex(1)
             -----------------------------------------
             \endverbatim
        */
        void formLexLSE(LexLSE& lexlse, Index& counter, Index ObjIndex)
        {
            Index CtrIndex;
            ConstraintType CtrType;

            if (ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP)
            {
                Index VarIndex;
                lexlse.setFixedVariablesCount(getActiveCtrCount());
                for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                {
                    CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                    VarIndex = getVarIndex(CtrIndex);             // CtrIndex       --> VarIndex
                    CtrType  = getActiveCtrType(CtrIndexActive);
                        
                    if (CtrType == LOWER_BOUND)
                        lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,0), LOWER_BOUND);
                    else if (CtrType == UPPER_BOUND)
                        lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,1), UPPER_BOUND);
                    else if (CtrType == EQUALITY_CONSTRAINT)
                        lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,1), EQUALITY_CONSTRAINT); // set equal to UPPER_BOUND by convention
                    else
                    {
                        throw Exception("CtrType is UNKNOWN");
                    }
                }
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE)
            {
                // FIXME: to implement ...
                throw Exception("WE SHOULD NOT BE HERE (YET)");
            }
            else if (ObjType == DEFAULT_OBJECTIVE)
            {
                RealScalar rhs;
                for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                {
                    CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                    CtrType  = getActiveCtrType(CtrIndexActive);

                    if (CtrType == EQUALITY_CONSTRAINT)
                    {
                        rhs = data.coeffRef(CtrIndex,nVar+1); // set equal to UPPER_BOUND by convention
                        lexlse.setCtrType(ObjIndex, CtrIndexActive, EQUALITY_CONSTRAINT);
                    } 
                    else if (CtrType == UPPER_BOUND)
                    {
                        rhs = data.coeffRef(CtrIndex,nVar+1); // set equal to UPPER_BOUND
                        lexlse.setCtrType(ObjIndex, CtrIndexActive, UPPER_BOUND);
                    }
                    else if (CtrType == LOWER_BOUND)
                    {
                        rhs = data.coeffRef(CtrIndex,nVar);   // set equal to LOWER_BOUND
                        lexlse.setCtrType(ObjIndex, CtrIndexActive, LOWER_BOUND);
                    }
                    else
                    {
                        throw Exception("UNKNOWN constraint type");
                    }                    

                    lexlse.setCtr(counter, data.row(CtrIndex).head(nVar), rhs);
                        
                    counter++;
                }
            }
        }

        /*
          \brief Form #dw
          
          \verbatim
          dw{inactive} =              0 - w{inactive}, where 0 corresponds to w{inactive}_star
          dw{active}   = w{active}_star - w{active}
          \endverbatim

          \return Squared norm of the step.
        */
        RealScalar formStep(dVectorType &dx)
        {
            RealScalar rhs;
            Index CtrIndex, lbIndex, ubIndex;
            ConstraintType CtrType;
            
            if (ObjType == DEFAULT_OBJECTIVE)
            {
                lbIndex = nVar;
                ubIndex = nVar+1;
             
                Adx = data.leftCols(nVar)*dx; // form Adx
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE || ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP)
            {
                lbIndex = 0;
                ubIndex = 1;
                
                for (Index k=0; k<nCtr; k++) // form Adx
                    Adx(k) = dx(VarIndex[k]); 
            }
            else
            {
                throw Exception("Unknown objective type");
            }
            
            // Only the first getActiveCtrCount() entries of wStar are initialized and used, so
            // wStar.setZero() is not necessary

            dw = -w;
            for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
            {
                CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                CtrType  = getActiveCtrType(CtrIndexActive);

                if (CtrType == EQUALITY_CONSTRAINT)
                    rhs = data.coeffRef(CtrIndex,ubIndex); // take upper bound by convention
                else if (CtrType == UPPER_BOUND)
                    rhs = data.coeffRef(CtrIndex,ubIndex);
                else if (CtrType == LOWER_BOUND)
                    rhs = data.coeffRef(CtrIndex,lbIndex);
                else
                    throw Exception("UNKNOWN constraint type"); // we should not be here
                
                wStar.coeffRef(CtrIndexActive) = Ax.coeffRef(CtrIndex) + Adx.coeffRef(CtrIndex) - rhs;

                dw.coeffRef(CtrIndex) += wStar.coeffRef(CtrIndexActive); // w{active}_star - w{active}
            }
            dwSquaredNorm = dw.squaredNorm();

            return dwSquaredNorm;
        }

        /**
           \brief Check for blocking constraints

           \param[out] CtrIndexBlocking Index of blocking constraint.
           \param[out] CtrTypeBlocking  Type of the blocking constraint.
           \param[out] alpha            scaling factor for the step.

           \return true if there are blocking constraints

           \verbatim
           -----------------------------------------
                [a', -1]*([x;w] + alpha*[dx;dw]) <= b : UPPER_BOUND
           b <= [a', -1]*([x;w] + alpha*[dx;dw])      : LOWER_BOUND
           -----------------------------------------
           den: a'*dx - dw
           num: b - a'*x + w
           ratio: num/den <-- this should be > 0 (in theory)
           
           UPPER_BOUND: check den > 0
           LOWER_BOUND: check den < 0

           alpha \in [0,1]
           \endverbatim
        */
        bool checkBlockingConstraints(Index &CtrIndexBlocking, ConstraintType &CtrTypeBlocking, RealScalar &alpha)
        {
            bool condition;
            RealScalar num, den, ratio, rhs, alpha_input = alpha;

            Index CtrIndex, lbIndex, ubIndex;
            ConstraintType CtrType;

            // I was using tol = 0, but if I have the same constraint appearing twice in one objective
            // (i.e., both RHS and LHS are the same) there is cycling. For tol = 1e-13 there is no cycling. 
            double tol = 1e-13; // FIXME: what to do with this tolerance?

            if (ObjType == DEFAULT_OBJECTIVE)
            {
                lbIndex = nVar;
                ubIndex = nVar+1;
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP || ObjType == SIMPLE_BOUNDS_OBJECTIVE)
            {
                lbIndex = 0;
                ubIndex = 1;
            }
            else
            {
                throw Exception("Unknown objective type");
            }
                        
            for (Index CtrIndexInactive=0; CtrIndexInactive<getInactiveCtrCount(); CtrIndexInactive++) // loop over inactive constraints
            {
                CtrIndex = getInactiveCtrIndex(CtrIndexInactive); // CtrIndexActive --> CtrIndex
                
                den = Adx.coeffRef(CtrIndex) - dw.coeffRef(CtrIndex);
                
                condition = false;
                if (den < -tol)     // LOWER_BOUND
                {
                    CtrType   = LOWER_BOUND;
                    rhs       = data.coeffRef(CtrIndex,lbIndex);
                    condition = true;
                }
                else if (den > tol) // UPPER_BOUND
                {
                    CtrType   = UPPER_BOUND;
                    rhs       = data.coeffRef(CtrIndex,ubIndex);
                    condition = true;
                }
                        
                if (condition)
                {
                    num   = rhs - Ax.coeffRef(CtrIndex) + w.coeffRef(CtrIndex);
                    ratio = num/den;

                    // ratio should always be positive (but just in case)
                    if (ratio < 0)
                        ratio = 0;
                            
                    // when den is very small, ratio will be grater than alpha (I don't expect numerical problems here)
                    if (ratio < alpha)
                    {
                        alpha = ratio;
                                
                        // store the current blocking constraint
                        CtrIndexBlocking = CtrIndex;
                        CtrTypeBlocking  = CtrType;
                    }
                }
                else
                {
                    //do nothing
                }
            }                
            
            return (alpha < alpha_input); // true if alpha has been modified
        }

        /**
           \brief Take a step alpha from #w in the direction #dw (and update #Ax)
           
           \param[in] alpha step scaling
        */
        void step(RealScalar alpha)
        {
            w  += alpha*dw;
            Ax += alpha*Adx;
        }
        
        // --------------------------------------------------------------------
        // set & get
        // --------------------------------------------------------------------
        
        /**
           \brief Returns #wStar
        */
        dVectorType& getOptimalResidual()
        {
            return wStar;
        }

        /**
           \brief Returns #w
        */
        dVectorType& getResidual()
        {
            return w;
        }

        /**
           \brief Returns #dw
        */
        dVectorType& get_dw()
        {
            return dw;
        }

        /**
           \brief Returns the number of active constraints
        */
        Index getActiveCtrCount() const
        {
            return WorkingSet.getActiveCtrCount();
        }

        /**
           \brief Returns the index of the k-th active constraint
        */
        Index getActiveCtrIndex(Index k) const
        {
            return WorkingSet.getActiveCtrIndex(k);
        }

        /**
           \brief Returns the type of the k-th active constraint
        */
        ConstraintType getActiveCtrType(Index k) const
        {
            return WorkingSet.getActiveCtrType(k);
        }

        /**
           \brief Returns the number of inactive constraints
        */
        Index getInactiveCtrCount() const
        {
            return WorkingSet.getInactiveCtrCount();
        }

        /**
           \brief Returns the index of the k-th inactive constraint
        */
        Index getInactiveCtrIndex(Index k) const
        {
            return WorkingSet.getInactiveCtrIndex(k);
        }

        /**
           \brief Returns the index storde in VarIndex[k] 
        */
        Index getVarIndex(Index k) const
        {
            return VarIndex(k);
        }

        /**
           \brief Returns a reference to VarIndex 
        */
        iVectorType& getVarIndex()
        {
            return VarIndex;
        }

        /**
           \brief Get objective type
        */
        ObjectiveType getObjType() const
        {
            return ObjType;
        }

        /**
           \brief Get number of constraints in objective
        */
        Index getDim() const
        {
            return nCtr;
        }

        /** 
            \brief Get objective data
        */
        MatrixType& getData()
        {
            return data;
        }

        /**
           \brief Returns true if the k-th constraint is active, otherwise returns false
        */                                        
        bool isActive(Index CtrIndex) const
        {
            return WorkingSet.isActive(CtrIndex);
        }

        /**
           \brief Set objective data (ObjType = DEFAULT_OBJECTIVE)
        */
        void setData(const MatrixType& data_)
        {
            data = data_;
        }

        /**
           \brief Set objective data + VarIndex (ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP or ObjType = SIMPLE_BOUNDS_OBJECTIVE)
        */
        void setData(Index *VarIndex_, const MatrixType& data_)
        {
            VarIndex = Eigen::Map<iVectorType>(VarIndex_,nCtr);
            data     = data_;
        }

        /**
           \brief Set w
        */
        void set_w(const dVectorType& w_)
        {
            w = w_;
            w_is_initialized = true;
        }

        /**
           \brief Modify upper or lower bound
        */
        void relax_bounds(Index CtrIndex, ConstraintType CtrType, RealScalar p)
        {
            Index lbIndex, ubIndex;

            // ----------------------------------------------------
            if (ObjType == DEFAULT_OBJECTIVE)
            {
                lbIndex = nVar;
                ubIndex = nVar+1;
            }
            else if (ObjType == SIMPLE_BOUNDS_OBJECTIVE_HP)
            {
                lbIndex = 0;
                ubIndex = 1;
            }
            else
            {
                throw Exception("Unknown objective type");
            }
            // ----------------------------------------------------
            if (CtrType == LOWER_BOUND)
            {
                data(CtrIndex,lbIndex) -= p; // relax lower-bound
            }
            else if (CtrType == UPPER_BOUND)
            {
                data(CtrIndex,ubIndex) += p; // relax upper-bound
            }
            else
            {
                throw Exception("Should not be here");
            }
            // ----------------------------------------------------
                
        }


        /**
           \brief Set objective data + VarIndex one by one (ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP or ObjType = SIMPLE_BOUNDS_OBJECTIVE)
        */
        void setData(Index k, Index VarIndex_, RealScalar bl, RealScalar bu)
        {
            VarIndex(k) = VarIndex_;
            data(k,0)   = bl;
            data(k,1)   = bu;
        }

        /**
           \brief Prints some fields

           \param[in] field description of field to print.
        */
        void print(const char * field) const
        {
            if (!strcmp(field, "WorkingSet"))
            {
                WorkingSet.print();
            }
            else if (!strcmp(field, "data"))
            {
                std::cout << "data = \n" << data << std::endl << std::endl;
                if (ObjType == SIMPLE_BOUNDS_OBJECTIVE)
                    std::cout << "VarIndex = {" << VarIndex.transpose() << "}"<< std::endl;
            }
            else if (!strcmp(field, "w"))
            {
                std::cout << "w = \n" << w << "\n dw  = \n" << dw << std::endl << std::endl;
            }
        }

    private:

        /**
           \brief Initializations
        */
        void initialize()
        {
            w.setZero();
            dw.setZero();

            dwSquaredNorm = 0;
        }

        // ------------------------------------------------------------------------------------------
        // fields
        // ------------------------------------------------------------------------------------------

        /**
           \brief Number of variables

           \note This field is a bit redundant as the number of variables is stored in LexLSI, but
           it is a bit inconvenient to always pass it as input argument.
        */
        Index nVar;

        /**
           \brief Number of constraints involved in objective
        */
        Index nCtr;

        /**
           \brief Objective type
        */
        ObjectiveType ObjType;

        /**
           \brief Vector of indexes of variables that are bounded (if ObjType =
           SIMPLE_BOUNDS_OBJECTIVE_HP or ObjType = SIMPLE_BOUNDS_OBJECTIVE)
        */
        iVectorType VarIndex;

        /**
           \brief Objective data
           
           \verbatim
           if ObjType = DEFAULT_OBJECTIVE         : data = [A, LowerBounds, UpperBounds]
           if ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP: data = [LowerBounds, UpperBounds], VarIndex is used
           if ObjType = SIMPLE_BOUNDS_OBJECTIVE   : data = [LowerBounds, UpperBounds], VarIndex is used
           \endverbatim
        */
        MatrixType data;

        /**
           \brief Working set
        */
        WorkingSetType WorkingSet;

        /**
           \brief Residual
        */
        dVectorType w;

        /**
           \brief Descent direction from #w
        */
        dVectorType dw;

        /**
           \brief When ObjType = DEFAULT_OBJECTIVE, stores A*x
        */
        dVectorType Ax;

        /**
           \brief When ObjType = DEFAULT_OBJECTIVE, stores A*dx
        */
        dVectorType Adx;

        /**
           \brief The first getActiveCtrCount() elements contain the optimal residual for the
           current getActiveCtrCount() equality constraints. The tail of wStar is never used.

           \note This variable is not used.
        */
        dVectorType wStar;

        /**
           \brief Squared norm of #dw
        */
        RealScalar dwSquaredNorm;

        /** 
            \brief If w_is_initialized == true, the function set_w(w) has been called.
        */
        bool w_is_initialized;
    };

} // END namespace LexLS

#endif // OBJECTIVE
