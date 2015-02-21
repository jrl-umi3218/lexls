#ifndef OBJECTIVE
#define OBJECTIVE

#include <utility.h>
#include <workingset.h>

namespace LexLS
{
    namespace internal
    {    
        /**
           \brief Defines an objective of a LexLSI problem
        */
        class Objective
        {
        public:

            /**
               \brief Default constructor
            */
            Objective(): 
                nCtr(0),
                obj_type(GENERAL_OBJECTIVE),
                regularization_factor(0.0),
                v0_is_specified(false){}

            /**
               \brief Resize the objective

               \param[in] nCtr_     Number of constraints in an objective.
               \param[in] nVar_     Number of variables (without the constraint violations).
               \param[in] obj_type_ Objective type
            */
            void resize(Index nCtr_, Index nVar_, ObjectiveType obj_type_)
            {
                obj_type = obj_type_;
                nCtr     = nCtr_;
                nVar     = nVar_;

                working_set.resize(nCtr);

                v.resize(nCtr);
                dv.resize(nCtr);

                Ax.resize(nCtr);
                Adx.resize(nCtr);

                if (obj_type == GENERAL_OBJECTIVE)
                {
                    data.resize(nCtr,nVar + 2); // [A,LowerBounds,UpperBounds]

                    lb_index = nVar; 
                    ub_index = nVar+1; 
                }
                else if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                {
                    data.resize(nCtr,2);         // [LowerBounds,UpperBounds]
                    var_index.resize(nCtr);

                    lb_index = 0; 
                    ub_index = 1; 
                }
                else
                {
                    throw Exception("Unknown objective type");
                }

                initialize();
            }

            void form_initial_working_set(dVectorType &x, 
                                          bool modify_x_guess_enabled,
                                          bool modify_type_active_enabled,
                                          bool modify_type_inactive_enabled)
            {
                Index CtrIndex;

                if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                {
                    //if (isActive(CtrIndex))
                    
                }
                else if (obj_type == GENERAL_OBJECTIVE)
                {
                    
                    
                }
            }

            /**
               \brief Initialize Ax
            */
            void initialize_Ax(dVectorType &x)
            {
                if (obj_type == GENERAL_OBJECTIVE)
                {
                    Ax = data.leftCols(nVar)*x;
                }
                else if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                {
                    for (Index k=0; k<nCtr; k++)
                    {
                        Ax(k) = x(var_index[k]); 
                    }
                }
            }

            /**
               \brief Given x, generate a v such that (x,v) is a feasible initial pair for the
               constraints involved in the objective 

               \note #Ax is initialized
            */
            void phase1(dVectorType &x)
            {
                initialize_Ax(x);

                Index CtrIndex;
                ConstraintActivationType CtrType;
            
                // put in function initialize_w0

                v  = Ax - 0.5*(data.col(lb_index)+data.col(ub_index)); 
            
                // overwrite v for the active constraints
                for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                {
                    CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                    CtrType  = getActiveCtrType(CtrIndexActive);
                
                    if (CtrType == CTR_ACTIVE_LB)
                        v.coeffRef(CtrIndex) = Ax.coeffRef(CtrIndex) - data.coeffRef(CtrIndex,lb_index);
                    else if (CtrType == CTR_ACTIVE_UB)
                        v.coeffRef(CtrIndex) = Ax.coeffRef(CtrIndex) - data.coeffRef(CtrIndex,ub_index);
                }
            
                // For inactive constraints, even when lb[i] <= Ax[i] <= ub[i], a nonzero v[i] can be
                // generated. So overwrite v[i] = 0.
                for (CtrIndex=0; CtrIndex<nCtr; CtrIndex++)
                {
                    if (!isActive(CtrIndex))
                    {
                        if (Ax.coeffRef(CtrIndex) >= data.coeffRef(CtrIndex,lb_index) && 
                            Ax.coeffRef(CtrIndex) <= data.coeffRef(CtrIndex,ub_index))
                        {
                            v.coeffRef(CtrIndex) = 0.0;
                        }
                    }
                }

            }

            /** 
                \brief Includes in the working set the constraint with index CtrIndex (and sets its type)

                \param[in] CtrIndex         Index of constraint in a given LexLSI objective
                \param[in] type             Type of the constraint to be included in the working set.
            */
            void activate(Index CtrIndex, ConstraintActivationType type)
            {
                if (CtrIndex >= nCtr)
                {
                    throw Exception("CtrIndex >= nCtr");
                }

                working_set.activate(CtrIndex, type);
            }

            /** 
                \brief Removes from the working set the constraint with index CtrIndexActive

                \param[in] CtrIndexActive Index of constraint in the working set in a given objective,
                i.e., Obj[ObjIndex].working_set.active[CtrIndexActive] will be removed.
            */
            void deactivate(Index CtrIndexActive)
            {
                if (CtrIndexActive >= getActiveCtrCount())
                {
                    throw Exception("CtrIndexActive >= number of active constraints");
                }

                working_set.deactivate(CtrIndexActive);
            }

            /** 
                \brief Form an LexLSE problem (using the current working set)

                \verbatim
                -----------------------------------------
                Example (handling of simple bounds)
                -----------------------------------------
                nVar = 6
                var_index = {3,5,4}
                active bounds = {2,0}
                inactive bounds = {1}
                type = {CTR_ACTIVE_LB, 
                CTR_ACTIVE_UB}
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
                ConstraintActivationType CtrType;

                if (obj_type == SIMPLE_BOUNDS_OBJECTIVE) // cannot be regularized
                {
                    Index VarIndex;
                    lexlse.setFixedVariablesCount(getActiveCtrCount());
                    for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                    {
                        CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                        VarIndex = getVarIndex(CtrIndex);             // CtrIndex       --> VarIndex
                        CtrType  = getActiveCtrType(CtrIndexActive);
                        
                        if (CtrType == CTR_ACTIVE_LB)
                        {
                            lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,0), CTR_ACTIVE_LB);
                        }
                        else if (CtrType == CTR_ACTIVE_UB)
                        {
                            lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,1), CTR_ACTIVE_UB);
                        }
                        else if (CtrType == CTR_ACTIVE_EQ)
                        {
                            lexlse.fixVariable(VarIndex, data.coeffRef(CtrIndex,1), CTR_ACTIVE_EQ); // set equal to upper bound by convention
                        }
                    }
                }
                else if (obj_type == GENERAL_OBJECTIVE)
                {
                    RealScalar rhs;
                    for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                    {
                        CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                        CtrType  = getActiveCtrType(CtrIndexActive);

                        if (CtrType == CTR_ACTIVE_EQ)
                        {
                            rhs = data.coeffRef(CtrIndex,nVar+1); // set equal to upper bound by convention
                            lexlse.setCtrType(ObjIndex, CtrIndexActive, CTR_ACTIVE_EQ);
                        } 
                        else if (CtrType == CTR_ACTIVE_UB)
                        {
                            rhs = data.coeffRef(CtrIndex,nVar+1); // set equal to upper bound
                            lexlse.setCtrType(ObjIndex, CtrIndexActive, CTR_ACTIVE_UB);
                        }
                        else if (CtrType == CTR_ACTIVE_LB)
                        {
                            rhs = data.coeffRef(CtrIndex,nVar);   // set equal to CTR_ACTIVE_LB
                            lexlse.setCtrType(ObjIndex, CtrIndexActive, CTR_ACTIVE_LB);
                        }

                        lexlse.setCtr(counter, data.row(CtrIndex).head(nVar), rhs);
                        
                        counter++;
                    }
                    lexlse.setRegularizationFactor(ObjIndex,regularization_factor);                
                }
            }           

            /*
              \brief Form #dv
          
              \verbatim
              dv{inactive} =              0 - v{inactive}, where 0 corresponds to v{inactive}_star
              dv{active}   = v{active}_star - v{active}
              \endverbatim
            */
            void formStep(dVectorType &dx)
            {
                RealScalar rhs;
                Index CtrIndex;
                ConstraintActivationType CtrType;

                if (obj_type == GENERAL_OBJECTIVE)
                {
                    Adx = data.leftCols(nVar)*dx; // form Adx
                }
                else if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                {
                    for (Index k=0; k<nCtr; k++) // form Adx
                    {
                        Adx(k) = dx(var_index[k]); 
                    }
                }

                dv = -v;
                for (Index CtrIndexActive=0; CtrIndexActive<getActiveCtrCount(); CtrIndexActive++)
                {
                    CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex
                    CtrType  = getActiveCtrType(CtrIndexActive);

                    if (CtrType == CTR_ACTIVE_EQ)
                    {
                        rhs = data.coeffRef(CtrIndex,ub_index); // take upper bound by convention
                    }
                    else if (CtrType == CTR_ACTIVE_UB)
                    {
                        rhs = data.coeffRef(CtrIndex,ub_index);
                    }
                    else if (CtrType == CTR_ACTIVE_LB)
                    {
                        rhs = data.coeffRef(CtrIndex,lb_index);
                    }
                
                    // v{active}_star - v{active}
                    dv.coeffRef(CtrIndex) += Ax.coeffRef(CtrIndex) + Adx.coeffRef(CtrIndex) - rhs;
                }
            }

            /**
               \brief Check for blocking constraints

               \param[out] CtrIndexBlocking Index of blocking constraint.
               \param[out] CtrTypeBlocking  Type of the blocking constraint.
               \param[out] alpha            scaling factor for the step.
               \param[in] tol_feasibility   used to determine whether a constraint has been violated.

               \return true if there are blocking constraints

               \verbatim
               -----------------------------------------
               [a', -1]*([x;v] + alpha*[dx;dv]) <= b : CTR_ACTIVE_UB
               b <= [a', -1]*([x;v] + alpha*[dx;dv])      : CTR_ACTIVE_LB
               -----------------------------------------
               den: a'*dx - dv
               num: b - a'*x + v
               ratio: num/den <-- this should be > 0 (in theory)
           
               CTR_ACTIVE_UB: check den > 0
               CTR_ACTIVE_LB: check den < 0

               alpha \in [0,1]
               \endverbatim
            */
            bool checkBlockingConstraints(Index &CtrIndexBlocking, 
                                          ConstraintActivationType &CtrTypeBlocking, 
                                          RealScalar &alpha,
                                          const RealScalar tol_feasibility)
            {
                bool condition;
                RealScalar num, den, ratio, rhs, alpha_input = alpha;

                Index CtrIndex;
                ConstraintActivationType CtrType;
                        
                for (Index CtrIndexInactive=0; CtrIndexInactive<getInactiveCtrCount(); CtrIndexInactive++) // loop over inactive constraints
                {
                    CtrIndex = getInactiveCtrIndex(CtrIndexInactive); // CtrIndexActive --> CtrIndex
                
                    den = Adx.coeffRef(CtrIndex) - dv.coeffRef(CtrIndex);
                
                    condition = false;
                    if (den < -tol_feasibility)     // CTR_ACTIVE_LB
                    {
                        CtrType   = CTR_ACTIVE_LB;
                        rhs       = data.coeffRef(CtrIndex,lb_index);
                        condition = true;
                    }
                    else if (den > tol_feasibility) // CTR_ACTIVE_UB
                    {
                        CtrType   = CTR_ACTIVE_UB;
                        rhs       = data.coeffRef(CtrIndex,ub_index);
                        condition = true;
                    }
                        
                    if (condition)
                    {
                        num   = rhs - Ax.coeffRef(CtrIndex) + v.coeffRef(CtrIndex);
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
               \brief Take a step alpha from #v in the direction #dv (and update #Ax)
           
               \param[in] alpha step scaling
            */
            void step(RealScalar alpha)
            {
                v  += alpha*dv;
                Ax += alpha*Adx;
            }
        
            // --------------------------------------------------------------------
            // set & get
            // --------------------------------------------------------------------

            /**
               \brief Returns #v
            */
            dVectorType& get_v()
            {
                return v;
            }

            /**
               \brief Returns #dv
            */
            dVectorType& get_dv()
            {
                return dv;
            }

            /**
               \brief Returns the number of active constraints
            */
            Index getActiveCtrCount() const
            {
                return working_set.getActiveCtrCount();
            }

            /**
               \brief Returns the index of the k-th active constraint
            */
            Index getActiveCtrIndex(Index k) const
            {
                return working_set.getActiveCtrIndex(k);
            }

            /**
               \brief Returns the type of the k-th active constraint
            */
            ConstraintActivationType getActiveCtrType(Index k) const
            {
                return working_set.getActiveCtrType(k);
            }

            /**
               \brief Returns the type of the k-th constraint
            */
            ConstraintActivationType getCtrType(Index k) const
            {
                return working_set.getCtrType(k);
            }

            /**
               \brief Returns the number of inactive constraints
            */
            Index getInactiveCtrCount() const
            {
                return working_set.getInactiveCtrCount();
            }

            /**
               \brief Returns the index of the k-th inactive constraint
            */
            Index getInactiveCtrIndex(Index k) const
            {
                return working_set.getInactiveCtrIndex(k);
            }

            /**
               \brief Returns the index storde in var_index[k] 
            */
            Index getVarIndex(Index k) const
            {
                return var_index(k);
            }

            /**
               \brief Get objective type
            */
            ObjectiveType getObjType() const
            {
                return obj_type;
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
            dMatrixType& getData()
            {
                return data;
            }

            /**
               \brief Returns true if the k-th constraint is active, otherwise returns false
            */
            bool isActive(Index CtrIndex) const
            {
                return working_set.isActive(CtrIndex);
            }
            
            /**
               \brief Set v0

               \note Use this function with caution (advanced initialization)
            */
            void set_v0(const dVectorType& v_)
            {
                v = v_;
                v0_is_specified = true;
            }

            /**
               \brief Modify upper or lower bound
            */
            void relax_bounds(Index CtrIndex, ConstraintActivationType CtrType, RealScalar p)
            {
                if (CtrType == CTR_ACTIVE_LB)
                {
                    data(CtrIndex,lb_index) -= p; // relax lower-bound
                }
                else if (CtrType == CTR_ACTIVE_UB)
                {
                    data(CtrIndex,ub_index) += p; // relax upper-bound
                }
                else
                {
                    throw Exception("Should not be here");
                }
            }

            /**
               \brief Set objective data (obj_type = GENERAL_OBJECTIVE)
            */
            void setData(const dMatrixType& data_)
            {
                data = data_;
            }

            /**
               \brief Set objective data + var_index (obj_type = SIMPLE_BOUNDS_OBJECTIVE)
            */
            void setData(Index *var_index_, const dMatrixType& data_)
            {
                var_index = Eigen::Map<iVectorType>(var_index_,nCtr);
                data      = data_;
            }

            /**
               \brief Set objective data + var_index one by one (obj_type = SIMPLE_BOUNDS_OBJECTIVE)
            */
            void setData(Index k, Index var_index_, RealScalar lb_, RealScalar ub_)
            {
                var_index(k) = var_index_;
                data(k,0)    = lb_;
                data(k,1)    = ub_;
            }

            /**
               \brief Set regularization factor
            */
            void setRegularization(RealScalar factor)
            {
                if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                {
                    printf("WARNING: setting a nonzero regularization factor has no effect on an objective of type SIMPLE_BOUNDS_OBJECTIVE. \n");
                }
                    
                regularization_factor = factor;
            }

            /**
               \brief Get the regularization factor
            */
            RealScalar getRegularization()
            {
                return regularization_factor;
            }

            /**
               \brief Prints some fields

               \param[in] field description of field to print.
            */
            void print(const char * field) const
            {
                if (!strcmp(field, "working_set"))
                {
                    working_set.print();
                }
                else if (!strcmp(field, "data"))
                {
                    std::cout << "data = \n" << data << std::endl << std::endl;
                    if (obj_type == SIMPLE_BOUNDS_OBJECTIVE)
                        std::cout << "var_index = {" << var_index.transpose() << "}"<< std::endl;
                }
                else if (!strcmp(field, "v"))
                {
                    std::cout << "v = \n" << v << "\n dv  = \n" << dv << std::endl << std::endl;
                }
            }

        private:

            /**
               \brief Initializations
            */
            void initialize()
            {
                v.setZero();
                dv.setZero();
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
               \brief Index of lower bound in the data

               \note It is different for general constraints and simple bounds
            */
            Index lb_index;

            /**
               \brief Index of upper bound in the data

               \note It is different for general constraints and simple bounds
            */
            Index ub_index;

            /**
               \brief Objective type
            */
            ObjectiveType obj_type;

            /**
               \brief Vector of indexes of variables that are bounded (if obj_type = SIMPLE_BOUNDS_OBJECTIVE)
            */
            iVectorType var_index;

            /**
               \brief Objective data
           
               \verbatim
               if obj_type = GENERAL_OBJECTIVE      : data = [A, LowerBounds, UpperBounds]
               if obj_type = SIMPLE_BOUNDS_OBJECTIVE: data = [LowerBounds, UpperBounds], var_index is used
               \endverbatim
            */
            dMatrixType data;

            /**
               \brief Working set
            */
            WorkingSet working_set;

            /**
               \brief Constraint violation
            */
            dVectorType v;

            /**
               \brief Descent direction from #v
            */
            dVectorType dv;

            /**
               \brief Stores A*x
            */
            dVectorType Ax;

            /**
               \brief Stores A*dx
            */
            dVectorType Adx;

            /*
              \brief Regularization factor for the current objective (default: 0.0)
            */
            RealScalar regularization_factor;

            /** 
                \brief If v0_is_initialized == true, the function set_v0(...) has been called.
            */
            bool v0_is_specified;
        };

    } // END namespace internal

} // END namespace LexLS

#endif // OBJECTIVE
