// Time-stamp: <2014-07-28 16:30:58 drdv>
#ifndef LEXLSI
#define LEXLSI

#include <lexlse.h>
#include <objective.h>
#include <cycling.h>

namespace LexLS
{
    /** 
        \brief Definition of a lexicographic least-squares problem with inequality constraints

        \todo When we solve a sequence of LexLSI problems we could specify the maximum size of the
        envisioned objectives so that we don't have to allocate memory online.

        \todo To use a structure containing the tolerances.
    */
    class LexLSI
    {
    public:

        // ---------------------------------------------------------------------------------
        // Constructors
        // ---------------------------------------------------------------------------------

        /** 
            \param[in] nVar_    Number of variables (only number of elements in x, and not in the residuals w)
            \param[in] nObj_    Number of objectives
            \param[in] ObjDim_  Number of constraints involved in each objective
            \param[in] ObjType_ Type of each objective
        */
        LexLSI(Index nVar_, Index nObj_, Index *ObjDim_, ObjectiveType *ObjType_):
            nVar(nVar_), 
            nObj(nObj_),
            iter(0),
            max_number_of_iterations(200),
            LinearDependenceTolerance(1e-12),
            DeactivationTolerance(1e-08),
            CyclingHandling(false),
            x0_is_initialized(false),
            status(TERMINATION_STATUS_UNKNOWN)
        {           
            resize(ObjDim_,ObjType_);
        }

        // ---------------------------------------------------------------------------------

        /** 
            \brief Adds a constraint to the working set (and sets its type)

            \param[in] ObjIndex         Index of objective.
            \param[in] CtrIndex         Index of constraint: Obj[ObjIndex].data.row(CtrIndex).
            \param[in] type             Type of the constraint.

            \note This function will be a part of the interface level and its purpose is to provide
            the initial working set.

            \todo Move all veification of inputs to the API level
        */
        void api_activate(Index ObjIndex, Index CtrIndex, ConstraintType type)
        {
            if (!Obj[ObjIndex].isActive(CtrIndex))
            {
                // which constraints are considered as EQUALITY_CONSTRAINT is determined internaly
                if (type == LOWER_BOUND || type == UPPER_BOUND)
                    activate(ObjIndex, CtrIndex, type, false);
                else // see setData(...)
                    std::cout << "WARNING: for the moment the user cannot define explicitly which constraints are of type EQUALITY_CONSTRAINT \n" << std::endl;
            }
        }
      
        /** 
            \brief Adds a constraint to the working set (and sets its type)

            \param[in] ObjIndex       Index of objective.
            \param[in] CtrIndex       Index of constraint: Obj[ObjIndex].data.row(CtrIndex).
            \param[in] type           Type of the constraint.
            \param[in] CountIteration if true, the iteration counter #iterAdd is incremented

            \note CountIteration = false is used when specifying the initial working set
        */
        void activate(Index ObjIndex, Index CtrIndex, ConstraintType type, bool CountIteration=true)
        {
            if (ObjIndex >= nObj)                
                throw Exception("ObjIndex >= nObj");

            Obj[ObjIndex].activate(CtrIndex, type);

            if (CountIteration)
                iterAdd++;
        }
        
        /** 
            \brief Removes a constraint from the working set

            \param[in] ObjIndex       Index of objective.
            \param[in] CtrIndexActive Index of constraint: Obj[ObjIndex].WorkingSet.active[CtrIndexActive].
        */
        void deactivate(Index ObjIndex, Index CtrIndexActive)
        {
            if (ObjIndex >= nObj)                
                throw Exception("ObjIndex >= nObj");

            Obj[ObjIndex].deactivate(CtrIndexActive);

            iterRemove++;
        }

        /**
           \brief Computes an initial feasible pair (x,w)
        */
        void phase1()
        {
            bool ActiveSet_is_initialized = false;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                if (Obj[ObjIndex].getActiveCtrCount() > 0)
                {
                    ActiveSet_is_initialized = true;
                    break;
                }
            }
         
            // --------------------------------------------------------
            // form x
            // --------------------------------------------------------
            if (ActiveSet_is_initialized)
            {                
                formLexLSE();                

                if (!x0_is_initialized)
                {
                    lexlse.factorize();
                    lexlse.solve();                    
                    x = lexlse.get_x();
                }
            }
            else
            {
                if (!x0_is_initialized)
                {
                    for (Index k=0; k<nVar; k++)
                        x(k) = 0.01; // set to something different from 0
                }
            }
            
            // --------------------------------------------------------
            // form w
            // --------------------------------------------------------
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                Obj[ObjIndex].phase1(x);

            // --------------------------------------------------------
            // form step (similar to formStep())
            // --------------------------------------------------------
            dx.setZero();
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                Obj[ObjIndex].formStep(dx);
            // --------------------------------------------------------

/*
            printf("x = [");
            for (Index k=0; k<nVar; k++)
                printf(" %f ",x[k]);
            printf("]\n\n");

            printf("dx = [");
            for (Index k=0; k<nVar; k++)
                printf(" %f ",dx[k]);
            printf("]\n\n");

            for (Index k=0; k<nObj; k++)
                Obj[k].print("w");
*/
        }
        
        /**
           \brief solve a LexLSI problem

           \return the termination reason
        */
        TerminationStatus solve()
        {
/*
            {
                std::ofstream pfile("./debug.dat", std::ios::out);
                pfile.close();
            }
*/
//            std::ofstream pfile("./debug.dat", std::ios::app);  
//            pfile.precision(15);

            phase1();

            while (1)
            {
                verifyWorkingSet();

//                pfile << "------------------------------------------------------" << "\n";
//                pfile << "iter = " << iter << "\n";
//                pfile << "------------------------------------------------------" << "\n";
//                pfile << getLambda() << "\n\n"; 

                //printf("iter = %d \n", iter);
                
                if ((status == PROBLEM_SOLVED) || (status == PROBLEM_SOLVED_CYCLING_HANDLING))
                {
                    //cycling_handler.print_counter();

                    break; // we are done ...
                }
                else
                {
                    if (getIterationsCount() > max_number_of_iterations)
                    {
                        status = MAX_NUMBER_OF_ITERATIONS_EXCEDED;
                        break;
                    }
                }
            }
//            pfile.close();

/*
            std::ofstream pfile("./debug.dat", std::ios::app);
            pfile << "solved (iter = " << iter << ")\n";
            pfile << "===================================================================" << "\n\n";
            pfile.close();
*/
            return status;
        }

        /**
           \brief Prints some fields
           
           \param[in] field description of field to print.

           \todo Remove this function.
        */
        void print(const char * field)
        {
            if (!strcmp(field, "WorkingSet"))
            {
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    Obj[ObjIndex].print("WorkingSet");
                std::cout << std::endl;
            }
            else if (!strcmp(field, "data"))
            {
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    std::cout << "--------------------------------------------------" << std::endl;
                    std::cout << "Objective[" << ObjIndex << "].";   
                    Obj[ObjIndex].print("data");
                }
                std::cout << std::endl;
            }
            else if (!strcmp(field, "iter"))
            {
                std::cout << "iter = "     << iter 
                          << " (ADD = "    << iterAdd 
                          << ", REMOVE = " << iterRemove 
                          << ", ACTIVE = " << getActiveCtrCount() << ")" << std::endl;
                std::cout << std::endl;
            }
            else if (!strcmp(field, "x"))
            {
                std::cout << "x = \n" << x << std::endl;
                std::cout << std::endl;
            }
            else if (!strcmp(field, "w"))
            {
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    std::cout << "w["<<ObjIndex<<"] = \n" << Obj[ObjIndex].getResidual() << std::endl;
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        }

        // ---------------------------------------------------------------------------------
        // set & get
        // ---------------------------------------------------------------------------------

        /**
           \brief Sets the initial value for the decision variable x
        */
        void set_x0(dVectorType &x0)
        {
            x = x0;
            x0_is_initialized = true;
        }

        /**
           \brief Sets the maximum number of iterations
        */
        void setMaxIter(Index max_number_of_iterations_)
        {
            max_number_of_iterations = max_number_of_iterations_;
        }

        /**
           \brief Sets tolerances
        */
        void setTolerance(RealScalar LinearDependenceTolerance_, 
                          RealScalar DeactivationTolerance_)
        {
            LinearDependenceTolerance = LinearDependenceTolerance_;
            DeactivationTolerance     = DeactivationTolerance_;

            lexlse.setTolerance(LinearDependenceTolerance);
        }

        /** 
            \brief Set data of objective ObjIndex (ObjType = DEFAULT_OBJECTIVE is assumed)

            \param[in] ObjIndex Index of objective
            \param[in] data     [A,LowerBounds,UpperBounds]
        */
        void setData(Index ObjIndex, const MatrixType& data)
        {
            if (ObjIndex >= nObj)                
                throw Exception("ObjIndex >= nObj");

            if (Obj[ObjIndex].getObjType() != DEFAULT_OBJECTIVE)
                throw Exception("ObjType = DEFAULT_OBJECTIVE is assumed");          

            if (Obj[ObjIndex].getDim() != data.rows())                
                throw Exception("Incorrect number of equations");

            // check bounds
            RealScalar bl, bu;
            for (Index CtrIndex=0; CtrIndex<Obj[ObjIndex].getDim(); CtrIndex++)
            {    
                bl = data.coeffRef(CtrIndex,nVar);
                bu = data.coeffRef(CtrIndex,nVar+1);
                             
                if (isEqual(bl,bu))
                    activate(ObjIndex,CtrIndex,EQUALITY_CONSTRAINT,false);
                else if (bl > bu)
                    throw Exception("(general) Lower bound is greater than upper bound.");   
            }

            Obj[ObjIndex].setData(data);
        }

        /** 
            \brief Set data of objective ObjIndex (ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP or ObjType = SIMPLE_BOUNDS_OBJECTIVE is assumed)

            \param[in] ObjIndex Index of objective
            \param[in] VarIndex Index variables subject to simple bounds
            \param[in] data     [LowerBounds,UpperBounds]
        */
        void setData(Index ObjIndex, Index *VarIndex, const MatrixType& data)
        {
            if (ObjIndex >= nObj)                
                throw Exception("ObjIndex >= nObj");

            if (Obj[ObjIndex].getObjType() != SIMPLE_BOUNDS_OBJECTIVE_HP && Obj[ObjIndex].getObjType() != SIMPLE_BOUNDS_OBJECTIVE)
                throw Exception("ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP or ObjType = SIMPLE_BOUNDS_OBJECTIVE is assumed");

            if (Obj[ObjIndex].getDim() != data.rows())                
                throw Exception("Incorrect number of equations");
            
            // check bounds
            RealScalar bl, bu;
            for (Index CtrIndex=0; CtrIndex<Obj[ObjIndex].getDim(); CtrIndex++)
            {    
                bl = data.coeffRef(CtrIndex,0);
                bu = data.coeffRef(CtrIndex,1);

                if (isEqual(bl,bu))
                    activate(ObjIndex,CtrIndex,EQUALITY_CONSTRAINT,false);
                else if (bl > bu)
                    throw Exception("(simple) Lower bound is greater than upper bound.");
            }

            // check whether VarIndex contains repeated indexes (VarIndex is not assumed to be sorted)
            for (Index k=0; k<Obj[ObjIndex].getDim(); k++)
                for (Index j=0; j<Obj[ObjIndex].getDim(); j++)
                    if ((VarIndex[k] == VarIndex[j]) && (j != k))
                        throw Exception("Elements of VarIndex are not unique.");   
            
            Obj[ObjIndex].setData(VarIndex, data);
        }

        /**
           \brief Sets #CyclingHandling
        */
        void setCyclingHandling(bool CyclingHandling_)
        {
            CyclingHandling = CyclingHandling_;
        }

        /** 
            \brief Return the (primal) solution vector
        */
        dVectorType& get_x()
        {
            return x;
        }

        dVectorType& getResidual(Index ObjIndex)
        {
            return Obj[ObjIndex].getResidual();
        }

        /** 
            \brief Outputs the Lagrange multipliers associated to the constraintes involved in all objectives

            \note The column corresponding to SIMPLE_BOUNDS_OBJECTIVE_HP is stored

            \note The order of the constraints in the active set is preserved. 

            \attention Note that L is returned by value.
        */
        MatrixType getLambda()
        {
            Index nActiveCtr = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                nActiveCtr += lexlse.getDim(ObjIndex); 
            
            MatrixType L = MatrixType::Zero(nActiveCtr,nObj);

            Index nMeaningful = lexlse.getFixedVariablesCount();
            for (Index ObjIndex=0; ObjIndex<nObj-ObjOffset; ObjIndex++) // Objectives of LexLSE
            {
                lexlse.ObjectiveSensitivity(ObjIndex);

                nMeaningful += lexlse.getDim(ObjIndex);
                L.col(ObjOffset + ObjIndex).head(nMeaningful) = lexlse.getWorkspace().head(nMeaningful);
            }

            return L;
        }

        /** 
            \brief Returns number of iterations in the active-set method
        */
        Index getIterationsCount() const
        {
            return iter;
        }

        /** 
            \brief Returns number of iterations during which a constraint has been added to the
            working set
        */
        Index getAddCount() const
        {
            return iterAdd;
        }

        /** 
            \brief Returns number of iterations during which a constraint has been removed to the
            working set
        */
        Index getRemoveCount() const
        {
            return iterRemove;
        }

        /** 
            \brief Returns number of active constraints
        */
        Index getActiveCtrCount() const
        {
            Index n = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                n += Obj[ObjIndex].getActiveCtrCount();

            return n;
        }

        /** 
            \brief Outputs the indexes and types of active constraints 
            (in the order they appear in the active set) for a given objective
        */
        void getActiveCtrOrder(Index ObjIndex, std::vector<Index>& indexes, std::vector<ConstraintType>& types) const
        {
            for (Index k=0; k<Obj[ObjIndex].getActiveCtrCount(); k++)
            {           
                indexes.push_back(Obj[ObjIndex].getActiveCtrIndex(k));
                types.push_back(Obj[ObjIndex].getActiveCtrType(k));
            }
        }

        /** 
            \brief Outputs the indexes of active constraints
            (in the order they appear in the active set) for a given objective
        */
        void getActiveCtrOrder(Index ObjIndex, std::vector<Index>& indexes) const
        {
            for (Index k=0; k<Obj[ObjIndex].getActiveCtrCount(); k++)
                indexes.push_back(Obj[ObjIndex].getActiveCtrIndex(k));
        }

        /** 
            \brief Outputs the indexes of active constraints for a given objective
        */
        void getActiveCtr(Index ObjIndex, std::vector<ConstraintType>& ctr_type) const
        {
            Index ind;
            Index dim = Obj[ObjIndex].getDim();
            ctr_type.resize(dim,CONSTRAINT_TYPE_UNKNOWN);
            for (Index k=0; k<Obj[ObjIndex].getActiveCtrCount(); k++)
            {   
                ind = Obj[ObjIndex].getActiveCtrIndex(k);
                ctr_type[ind] = Obj[ObjIndex].getActiveCtrType(k);
            }
        }

        /** 
            \brief Returns number of objectives
        */
        Index getObjectivesCount() const
        {
            return nObj;
        }

        /** 
            \brief Returns number of constraints in objective ObjIndex

            \param[in] ObjIndex Index of objective
        */
        Index getObjDim(Index ObjIndex) const
        {
            return Obj[ObjIndex].getDim();
        }
        
        /** 
            \brief Reset the LexLSI problem (in order to resolve it)
        */
        void reset()
        {
            Index ActiveCtrCount, counter;

            // Deactivate all active constraints that are not of type EQUALITY_CONSTRAINT.
            // Note that WorkingSet.inactive would not be sorted and the order of constraints would, in general, change
            // when the same problem is resolved multiple times (using reset()) - this should not cause a problem.  
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {         
                counter = 0;
                ActiveCtrCount = Obj[ObjIndex].getActiveCtrCount();
                for (Index CtrIndexActive=0; CtrIndexActive<ActiveCtrCount; CtrIndexActive++)
                {
                    if (Obj[ObjIndex].getActiveCtrType(counter) != EQUALITY_CONSTRAINT)
                    {
                        deactivate(ObjIndex, counter);
                    }
                    else
                    {
                        counter++;
                    }
                }
            }
            
            initialize();

            iter   = 0;
            status = TERMINATION_STATUS_UNKNOWN;
        }
        
    private:
        
        /** 
            \brief Resize LexLSI problem

            \param[in] ObjDim_  Number of constraints involved in each objective
            \param[in] ObjType_ Type of each objective
        */
        void resize(Index *ObjDim_, ObjectiveType *ObjType_)
        {
            ObjOffset = 0;
            if (ObjType_[0] == SIMPLE_BOUNDS_OBJECTIVE_HP) // If only simple bounds in first objective of LexLSI
                ObjOffset = 1;

            // In LexLSE, fixed variables are handled separately and are not defined as an objective
            // ObjDim_ + ObjOffset is pointer arithmetic
            lexlse.resize(nVar, nObj - ObjOffset, ObjDim_ + ObjOffset);

            nActive.resize(nObj);
            Obj.resize(nObj); 
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                Obj[ObjIndex].resize(ObjDim_[ObjIndex],nVar,ObjType_[ObjIndex]);

            x.resize(nVar);
            dx.resize(nVar); 

            initialize();
        }

        /** 
            \brief Initializations
        */
        void initialize()
        {            
            iterAdd    = 0;
            iterRemove = 0;

            x.setZero();
            dx.setZero();   
        }

        /** 
            \brief Form an LexLSE problem (using the current working set)
        */
        void formLexLSE()
        {
            // ObjInfo.FirstRowIndex has to be initialized before I start setting CtrType in formLexLSE below
            Index counter = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                nActive(ObjIndex) = Obj[ObjIndex].getActiveCtrCount();
            lexlse.setObjDim(&nActive(0)+ObjOffset);
            
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                Obj[ObjIndex].formLexLSE(lexlse, counter, ObjIndex-ObjOffset);
        }

        /**
           \brief Form the step (dx,dw) from the current iterate and compute the step length StepLength

           \return the squared norm of (dx,dw)
        */
        RealScalar formStep()
        {            
            dx = lexlse.get_x() - x;

            RealScalar StepLength = dx.squaredNorm();
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                StepLength += Obj[ObjIndex].formStep(dx);
            
            return StepLength;
        }

        /**
           \brief Check for blocking constraints

           \param[out] ObjIndexBlocking Index of objective with the blocking constraint (if such exists).
           \param[out] CtrIndexBlocking Index of blocking constraint (in LexLSI objective ObjIndexBlocking).
           \param[out] CtrTypeBlocking  Type of the blocking constraint.
           \param[out] alpha            scaling factor for the step.

           \return true if there are blocking constraints
        */      
        bool checkBlockingConstraints(Index &ObjIndexBlocking, 
                                      Index &CtrIndexBlocking, 
                                      ConstraintType &CtrTypeBlocking, 
                                      RealScalar &alpha)
        {
            alpha = 1;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                if (Obj[ObjIndex].checkBlockingConstraints(CtrIndexBlocking, CtrTypeBlocking, alpha))
                    ObjIndexBlocking = ObjIndex;
            
            if (alpha < 1)
                return true; // there are blocking constraints
            else
                return false;
        }

        /**
           \brief Finds active constraints that should be removed from the working set

           \param[out] ObjIndex2Remove Index of objective from which to remove a constraint
           \param[out] CtrIndex2Remove Index of constraint in the working set of objective ObjIndex2Remove to remove

           \return true if there are constraints to remove
        */
        bool findActiveCtr2Remove(Index &ObjIndex2Remove, Index &CtrIndex2Remove)
        {
            bool DescentDirectionExists = false;
            RealScalar tol, residual_norm;
            Index dim;
            
            for (Index ObjIndex=0; ObjIndex<nObj-ObjOffset; ObjIndex++) // loop over objectives of LexLSE problem
            {
                // --------------------------------------------------
                // determine exit tolerance for Lagrange multipliers
                // --------------------------------------------------
                if (0)
                {
                    dim = 0;
                    residual_norm = 0;
                    for (Index k=ObjIndex; k<nObj-ObjOffset; k++) 
                    {
                        dim           += Obj[k].getDim();
                        residual_norm += Obj[k].get_wStarSquaredNorm();
                    }
                    tol = dim*DeactivationTolerance;
                    if (residual_norm > 1)
                        tol *= residual_norm;
                }
                else
                {
                    tol = DeactivationTolerance;
                }
                // --------------------------------------------------
                
                DescentDirectionExists = lexlse.ObjectiveSensitivity(ObjIndex, CtrIndex2Remove, ObjIndex2Remove, tol);
                if (DescentDirectionExists)
                    break;
            }
            
            // Note that when the first objective of LexLSI is of type SIMPLE_BOUNDS_OBJECTIVE_HP,
            // and if a constraint is to be removed from it, ObjIndex2Remove = -1 (see end of LexLSE.ObjectiveSensitivity(...)).
            ObjIndex2Remove += ObjOffset; // objective of LexLSI problem

            return DescentDirectionExists;
        }

        /**
           \brief One iteration of an active-set method
        */        
        void verifyWorkingSet()
        {
            // ----------------------------------------------------------------------
            Index ObjIndex2Manipulate, CtrIndex2Manipulate;
            ConstraintType CtrType2Manipulate = CONSTRAINT_TYPE_UNKNOWN;

            bool normalIteration = true;         
            OperationType operation = UNDEFINED;
            ConstraintIdentifierType ConstraintIdentifier;

            RealScalar alpha;
            // ----------------------------------------------------------------------

/*
            printf("----------------------------------------------\n\n");
            printf("iter = %d \n", iter);
            printf("----------------------------------------------\n\n");

            printf("x = \n");
            for (Index i=0; i<nVar; i++)
            {
                printf("%f \n", x[i]);
            }
            printf("\n");

            printf("dx = \n");
            for (Index i=0; i<nVar; i++)
            {
                printf("%f \n", dx[i]);
            }
            printf("\n");

            Obj[0].print("w");
*/

            if (iter != 0) // iter == 0 is handled in phase1()
            {
                formLexLSE();
                
                lexlse.factorize();
                lexlse.solve();

                formStep();
            }
            else // if iter == 0
            {
                if (x0_is_initialized)
                {
                    normalIteration = false;

                    /// @todo Commented out.
                    //printf("\n\nWARNING: normalIteration = false \n\n");
                }
            }

            if (checkBlockingConstraints(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate, alpha))
            {
/*
                printf("(%4d)   ADD(%d): obj = %d, ctr = %3d (dim = %3d, alpha = %+e)\n", 
                       iter, 
                       CtrType2Manipulate, 
                       ObjIndex2Manipulate, 
                       CtrIndex2Manipulate, 
                       Obj[ObjIndex2Manipulate].getActiveCtrCount(), 
                       alpha);
*/
                if (CyclingHandling)
                {
                    ConstraintIdentifier.set(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate);
                }

                operation = ADD_CONSTRAINT;
                activate(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate);
            }
            else
            {
                if (normalIteration) 
                {
                    if (findActiveCtr2Remove(ObjIndex2Manipulate, CtrIndex2Manipulate))
                    {
/*
                        printf("(%4d)REMOVE(%d): obj = %d, ctr = %3d (dim = %3d)\n", 
                               iter, 
                               Obj[ObjIndex2Manipulate].getActiveCtrType(CtrIndex2Manipulate),
                               ObjIndex2Manipulate, 
                               Obj[ObjIndex2Manipulate].getActiveCtrIndex(CtrIndex2Manipulate), 
                               Obj[ObjIndex2Manipulate].getActiveCtrCount()-1); // the constraint is removed below
*/
                        if (CyclingHandling)
                        {
                            ConstraintIdentifier.set(ObjIndex2Manipulate, 
                                                     Obj[ObjIndex2Manipulate].getActiveCtrIndex(CtrIndex2Manipulate), 
                                                     Obj[ObjIndex2Manipulate].getActiveCtrType(CtrIndex2Manipulate));
                        }

                        operation = REMOVE_CONSTRAINT;
                        deactivate(ObjIndex2Manipulate, CtrIndex2Manipulate); 
                    }
                    else
                    {
                        status = PROBLEM_SOLVED;
                    }
                }
            }

            //printf("alpha = %f, status = %d\n", alpha, status);
            
            if (alpha > 0) // take a step
            {
                x += alpha*dx;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    Obj[ObjIndex].step(alpha);
            }

            if (CyclingHandling && operation != UNDEFINED)
                status = cycling_handler.update(operation, ConstraintIdentifier, Obj, iter, false);

            iter++;
        }

        // ==================================================================
        // definition of scalars
        // ==================================================================

        /** 
            \brief Number of decision variables #x

            \note If we consider the problem: minimize_{x,w} norm(w,2)^2, subject to A*x - b = w,
            then clearly w is a decision variable as well, but we could always think of this propblem in
            terms of: minimize_{x} norm(A*x-b,2)^2.
        */
        Index nVar;  
    
        /** 
            \brief Number of objectives
        */
        Index nObj;

        /** 
            \brief When the objective with highest priority of LexLSI has only simple bounds (i.e.,
         its ObjType = SIMPLE_BOUNDS_OBJECTIVE_HP), the number of objectives in LexLSI and LexLSE differ
         with 1 because fixed variables are not treated as an objective in LexLSE.
        */
        Index ObjOffset;

        /*
          \brief Number of iterations
        */
        Index iter;

        /*
          \brief Number of iterations during which a constraint was added
        */
        Index iterAdd;

        /*
          \brief Number of iterations during which a constraint was removed
        */
        Index iterRemove;

        /*
          \brief Maximum number of iterations
        */
        Index max_number_of_iterations;

        /** 
            \brief Linear dependence tolerance (used when solving an LexLSE problem) 
        */
        RealScalar LinearDependenceTolerance;

        /** 
            \brief Sensitivity tolerance (used when determining the lexicographic sign of Lagrange multipliers)
        */
        RealScalar DeactivationTolerance;

        /** 
            \brief If CyclingHandling == true, cycling handling is performed
        */
        bool CyclingHandling;

        /** 
            \brief If x0_is_initialized == true, the function set_x0(dVectorType &x0) has been
            called and x0 has been initialized.

            \note This is later used in phase1().
        */
        bool x0_is_initialized;
    
        // ==================================================================
        // definition of vectors
        // ==================================================================

        /** 
            \brief The current value of the decision variables - not including the residual
        */
        dVectorType x;

        /** 
            \brief The current descent direction from #x
        */
        dVectorType dx;

        /** 
            \brief Number of active constraints in each objective
            
            \note This variable is used for convenience. 
        */
        iVectorType nActive;
        
        // ==================================================================
        // other definitions
        // ==================================================================

        /** 
            \brief Provides information about the reson for termination
        */
        TerminationStatus status;

        /** 
            \brief Handles the lexicographic least-squares problem with equality constraints

            \note This instance of LexLSE is used to solve multiplie problems - it is initialized
            with the largest expected problem dimensions. 
        */       
        LexLSE lexlse;

        /** 
            \brief Vector of objectives
        */       
        std::vector<Objective> Obj;

        /** 
            \brief Handle cycling
        */       
        CyclingHandlerType cycling_handler;

        /** 
            \brief Vector used to store constraints that have been added to the working set since
            the last time the step alpha = 0.
        */       
        //std::vector<ConstraintIdentifierType> CyclingArray;

        /** 
            \brief The same as CyclingArray but for the previous series of constraint additions.

            \note Essentially CyclingArrayPrevious stores CyclingArray before it is cleared
        */       
        //std::vector<ConstraintIdentifierType> CyclingArrayPrevious;

    }; // END class LexLSI 

} // END namespace LexLS 

#endif // LEXLSE
