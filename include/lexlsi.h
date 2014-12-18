// Time-stamp: <2014-12-18 17:51:22 drdv>
#ifndef LEXLSI
#define LEXLSI

#include <lexlse.h>
#include <objective.h>
#include <cycling.h>

namespace LexLS
{
    class LexLSIParameters
    {
    public:
        /**
           \brief Maximum number of iterations
        */
        Index max_number_of_iterations;

        /** 
            \brief Tolerance: linear dependence (used when solving an LexLSE problem)
        */
        RealScalar tolLinearDependence;

        /** 
            \brief Tolerance: absolute value of Lagrange multiplier to be considered with "wrong" sign
        */
        RealScalar tolWrongSignLambda;

        /** 
            \brief Tolerance: absolute value of Lagrange multiplier to be considered with "correct" sign
        */
        RealScalar tolCorrectSignLambda;

        /** 
            \brief If CyclingHandling == true, cycling handling is performed
        */
        bool CyclingHandling;
 
        RegularizationType regularizationType;

        Index regularizationMaxIterCG;

        Index cycling_max_counter;
        double cycling_relax_step;

        // use the real residual when computing sensitivity
        bool realSensitivityResidual;

        std::string output_file_name; // drdv: do I need to initialize?

        LexLSIParameters()
        {
            setDefaults();
        }

        void setDefaults()
        {
            max_number_of_iterations = 200;

            tolLinearDependence     = 1e-12;
            tolWrongSignLambda      = 1e-08;
            tolCorrectSignLambda    = 1e-12;

            CyclingHandling         = false;
            cycling_max_counter     = 50;
            cycling_relax_step      = 1e-08;

            regularizationType      = REGULARIZATION_NONE;

            regularizationMaxIterCG = 10;

            realSensitivityResidual = false;
        }
    };



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
            x0_is_initialized(false),
            status(TERMINATION_STATUS_UNKNOWN)
        {           
            parameters.setDefaults();
            setParameters(parameters);

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
            bool active_constraints_exist = false;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                // we would enter even if there are only equality constraints
                if (Obj[ObjIndex].getActiveCtrCount() > 0)
                {
                    active_constraints_exist = true;
                    break;
                }
            }
         
            // --------------------------------------------------------
            // form x
            // --------------------------------------------------------
            if (active_constraints_exist)
            {
                formLexLSE();                
                
                if (!x0_is_initialized)
                {
                    lexlse.factorize();
                    lexlse.solve();                    
                    x = lexlse.get_x();

                    numberOfFactorizations++;
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
        }
        
        /**
           \brief solve a LexLSI problem

           \return the termination reason
        */
        TerminationStatus solve()
        {
            OperationType operation;

            phase1();

            if (!parameters.output_file_name.empty())
                outputStuff(parameters.output_file_name.c_str(), UNDEFINED, true);

            while (1)
            {
                operation = verifyWorkingSet();

                if (!parameters.output_file_name.empty())
                    outputStuff(parameters.output_file_name.c_str(), operation);
              
                if ((status == PROBLEM_SOLVED) || (status == PROBLEM_SOLVED_CYCLING_HANDLING))
                {
                    break; // we are done ...
                }
                else
                {
                    if (numberOfFactorizations => parameters.max_number_of_iterations)
                    {
                        status = MAX_NUMBER_OF_ITERATIONS_EXCEDED;
                        break;
                    }
                }
            }
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
           \brief Sets the residual for objective k
        */
        void set_w(Index ObjIndex, dVectorType &w)
        {
            Obj[ObjIndex].set_w(w);
        }

        /**
           \brief Sets parameters
        */
        void setParameters(const LexLSIParameters &parameters_)
        {
            parameters = parameters_;

            lexlse.setTolerance(parameters.tolLinearDependence);
            lexlse.setRegularizationType(parameters.regularizationType);
            lexlse.setRegularizationMaxIterCG(parameters.regularizationMaxIterCG);
            lexlse.setRealSensitivityResidual(parameters.realSensitivityResidual);

            if (parameters.CyclingHandling)
            {
                cycling_handler.set_max_counter(parameters.cycling_max_counter);
                cycling_handler.set_relax_step(parameters.cycling_relax_step);
            }
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
            \brief Set (a non-negative) regularization factor for objective ObjIndex

            \note Regularization of an objective of type SIMPLE_BOUNDS_OBJECTIVE_HP is not performed
        */        
        void setRegularization(Index ObjIndex, RealScalar RegularizationFactor)
        {
            // @todo: check whether ObjIndex and RegularizationFactor make sense. 

            regularization(ObjIndex) = RegularizationFactor;
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
        Index getNumberOfFactorizations() const
        {
            return numberOfFactorizations;
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
        /// \brief Parameters of the solver.
        LexLSIParameters parameters;


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

            regularization.resize(nObj);

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
            iter       = 0;
            iterAdd    = 0;
            iterRemove = 0;
            numberOfFactorizations = 0;

            x.setZero();
            dx.setZero();

            regularization.setZero(); // by default there is no regularization
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
            
            // skip the first regularization factor if the first objective is of type SIMPLE_BOUNDS_OBJECTIVE_HP
            for (Index ObjIndex=0; ObjIndex<nObj-ObjOffset; ObjIndex++)
                lexlse.setRegularization(ObjIndex,regularization(ObjIndex + ObjOffset));

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
            for (Index ObjIndex=0; ObjIndex<nObj-ObjOffset; ObjIndex++) // loop over objectives of LexLSE problem
            {
                // The real residual Obj[ObjIndex].getOptimalResidual() is an input but it might not be used inside the function. 
                // In fact it is better to use the residual computed from the factorization in lexlse because 
                // I have observed less cycling. I am just testing some stuff with the the real residual 
                // (don't use it if you don't know what you are doing!)
                DescentDirectionExists = lexlse.ObjectiveSensitivity(ObjIndex, 
                                                                     CtrIndex2Remove, 
                                                                     ObjIndex2Remove, 
                                                                     parameters.tolWrongSignLambda,
                                                                     parameters.tolCorrectSignLambda,
                                                                     Obj[ObjIndex].getOptimalResidual());

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
        OperationType verifyWorkingSet()
        {
            // ----------------------------------------------------------------------
            Index ObjIndex2Manipulate, CtrIndex2Manipulate;
            ConstraintType CtrType2Manipulate = CONSTRAINT_TYPE_UNKNOWN;

            bool normalIteration = true;         
            OperationType operation = UNDEFINED;
            ConstraintIdentifierType ConstraintIdentifier;

            RealScalar alpha;
            // ----------------------------------------------------------------------

            if (iter != 0) // iter == 0 is handled in phase1()
            {
                formLexLSE();
                
                lexlse.factorize();
                lexlse.solve();

                formStep();

                numberOfFactorizations++;
            }
            else // if iter == 0
            {
                if (x0_is_initialized)
                {
                    normalIteration = false;
                }
            }

            if (checkBlockingConstraints(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate, alpha))
            {
                if (parameters.CyclingHandling)
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
                        if (parameters.CyclingHandling)
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

            if (alpha > 0) // take a step
            {
                x += alpha*dx;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    Obj[ObjIndex].step(alpha);
            }

            if (parameters.CyclingHandling && operation != UNDEFINED)
                status = cycling_handler.update(operation, ConstraintIdentifier, Obj, iter, false);

            iter++;

            return operation;
        }

        /**
           \brief Outputs resiadual norm to file
        */
        void outputStuff(const char *file_name, OperationType operation, bool flag_clear_file = false)
        {
            // clear the content of the file
            if (flag_clear_file)
            {
                std::ofstream file(file_name, std::ios::out | std::ofstream::trunc);
                file.close();
            }

            std::ofstream file(file_name, std::ios::out | std::ios::app);
            file.precision(15);

            file << "operation_("<<iter+1<<") = " << operation << ";\n"; 

            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                dVectorType w_ = Obj[ObjIndex].getResidual();

                file << "w_{"<<ObjIndex+1<<"}(:,"<<iter+1<<") = [ "; 
                for (Index k=0; k<Obj[ObjIndex].getDim(); k++)
                {
                    file << w_(k) << " ";
                }
                file << "]';\n";
            }

            file << "nw_(:,"<<iter+1<<") = [ "; 
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                file << std::sqrt(Obj[ObjIndex].getResidualSquaredNorm()) << " "; 
            file << "]';\n";

            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                file << "a_{"<<ObjIndex+1<<"}(:,"<<iter+1<<") = [ "; 
                for (Index k=0; k<Obj[ObjIndex].getDim(); k++)
                    file << (Index) Obj[ObjIndex].getCtrType(k) << " ";
                file << "]';\n";
            }

            file << "x_(:,"<<iter+1<<") = [ "; 
            for (Index k=0; k<nVar; k++)
                file << x(k) << " "; 
            file << "]'; \n"; 

            file << "\n";

            file.close();
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
          \brief Number of factorization
        */
        Index numberOfFactorizations;

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
            \brief A heuristic regularization for each objective

            \note The number of elements is equal to the number of objectives (having a non-zero
            regularization for objective of type SIMPLE_BOUNDS_OBJECTIVE_HP has no effect as it is
            ignored)
        */
        dVectorType regularization; 

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

    }; // END class LexLSI 

} // END namespace LexLS 

#endif // LEXLSE
