#ifndef LEXLSI
#define LEXLSI

#include <lexlse.h>
#include <objective.h>
#include <cycling.h>

namespace LexLS
{
    namespace internal
    {
        /**
            \brief Definition of a lexicographic least-squares problem with inequality constraints

            \todo When we solve a sequence of LexLSI problems we could specify the maximum size of the
            envisioned objectives so that we don't have to allocate memory online.

            \todo When the first objective is of type SIMPLE_BOUNDS_OBJECTIVE, the user has to
            specify the indexes of the variables and their corresponding bounds. To monitor whether
            the bounds for the i-th variables have already been specified and if the user sets them
            again to give an error message.

            \todo When solving a sequence of inequality constrained problem we could reuse some of
            the memory.

            \todo When we solve an equality constrained problem with lexlsi, the working_set_log is of size 0x1. This
            is because we don't push_back in it (which is normal since there are neither blocking constraints,
            nor Lagrange multipliers to drop). But what if I am interested only in the rank field! Maybe we
            should separate rank?

            \todo write a function to deactivate weakly active inequality constraints (when requested
            by the user). Reason: hot starting, maybe we can write equate oc QPs by discriminating
            equality and inequality constraints.
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
                x_guess_is_specified(false),
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
                \param[in] CtrIndex         Index of constraint: objectives[ObjIndex].data.row(CtrIndex).
                \param[in] type             Type of the constraint.

                \note This function will be a part of the interface level and its purpose is to provide
                the initial working set.

                \todo Move all veification of inputs to the API level
            */
            void api_activate(Index ObjIndex, Index CtrIndex, ConstraintActivationType type)
            {
                if (!objectives[ObjIndex].isActive(CtrIndex))
                {
                    // which constraints are considered as CTR_ACTIVE_EQ is determined internaly
                    if (type == CTR_ACTIVE_LB || type == CTR_ACTIVE_UB)
                    {
                        activate(ObjIndex, CtrIndex, type, false);
                    }
                    else // see setData(...)
                    {
                        std::cout << "WARNING: the user cannot define explicitly which constraints are of type CTR_ACTIVE_EQ \n" << std::endl;
                    }
                }
            }

            /**
                \brief Adds a constraint to the working set (and sets its type)

                \param[in] ObjIndex        Index of objective.
                \param[in] CtrIndex        Index of constraint: objectives[ObjIndex].data.row(CtrIndex).
                \param[in] type            Type of the constraint.
                \param[in] CountActivation if true, the iteration counter #nActivations is incremented

                \note CountActivation = false is used when specifying the initial working set
            */
            void activate(Index ObjIndex, Index CtrIndex, ConstraintActivationType type, bool CountActivation=true)
            {
                if (ObjIndex >= nObj)
                {
                    throw Exception("ObjIndex >= nObj");
                }

                WS.push_back(ConstraintInfo(ObjIndex,CtrIndex));

                objectives[ObjIndex].activate(CtrIndex, type);

                if (CountActivation)
                {
                    nActivations++;

                    if (objectives[ObjIndex].isZeroNormal(CtrIndex))
                    {
                        printf("WARNING: activated inequality constraint (0*x = b): (obj_index = %d, ctr_index = %d) \n", ObjIndex, CtrIndex);
                    }
                }
            }

            /**
                \brief Removes a constraint from the working set

                \param[in] ObjIndex          Index of objective.
                \param[in] CtrIndexActive    Index of constraint: objectives[ObjIndex].working_set.active[CtrIndexActive].
            */
            void deactivate(Index ObjIndex, Index CtrIndexActive)
            {
                if (ObjIndex >= nObj)
                {
                    throw Exception("ObjIndex >= nObj");
                }

                // -----------------------------------------------------------
                std::vector<ConstraintInfo>::iterator it;
                it = std::find(WS.begin(),
                               WS.end(),
                               ConstraintInfo(ObjIndex,
                                              objectives[ObjIndex].getActiveCtrIndex(CtrIndexActive)));
                WS.erase(it);
                // -----------------------------------------------------------

                objectives[ObjIndex].deactivate(CtrIndexActive);

                nDeactivations++;
            }

            /**
               \brief solve a LexLSI problem

               \return the termination reason
            */
            TerminationStatus solve()
            {
                OperationType operation;

                if (parameters.use_phase1_v0)
                {
                    phase1_v0();
                }
                else
                {
                    phase1();
                }

                if (!parameters.output_file_name.empty())
                {
                    outputStuff(parameters.output_file_name.c_str(), OPERATION_UNDEFINED, true);
                }

                while (1)
                {
                    operation = verifyWorkingSet();

                    if (!parameters.output_file_name.empty())
                    {
                        outputStuff(parameters.output_file_name.c_str(), operation);
                    }

                    if ((status == PROBLEM_SOLVED) || (status == PROBLEM_SOLVED_CYCLING_HANDLING))
                    {
                        break; // we are done ...
                    }
                    else
                    {
                        if (nFactorizations >= parameters.max_number_of_factorizations)
                        {
                            status = MAX_NUMBER_OF_FACTORIZATIONS_EXCEDED;
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
                if (!strcmp(field, "working_set"))
                {
                    for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    {
                        objectives[ObjIndex].print("working_set");
                    }
                    std::cout << std::endl;
                }
                else if (!strcmp(field, "data"))
                {
                    for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    {
                        std::cout << "--------------------------------------------------" << std::endl;
                        std::cout << "Objectives[" << ObjIndex << "].";
                        objectives[ObjIndex].print("data");
                    }
                    std::cout << std::endl;
                }
                else if (!strcmp(field, "nIterations"))
                {
                    std::cout << "nIterations = " << nIterations
                              << " (ADD = "       << nActivations
                              << ", REMOVE = "    << nDeactivations
                              << ", FACTOR = "    << nFactorizations
                              << ", ACTIVE = "    << getActiveCtrCount() << ")" << std::endl;
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
                        std::cout << "w["<<ObjIndex<<"] = \n" << objectives[ObjIndex].get_v() << std::endl;
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
                x_guess_is_specified = true;
            }

            /**
               \brief Sets the residual for objective k

               \todo Check the validity of the hot-start
            */
            void set_v0(Index ObjIndex, dVectorType &v0_)
            {
                objectives[ObjIndex].set_v0(v0_);
            }

            /**
               \brief Sets parameters
            */
            void setParameters(const ParametersLexLSI &parameters_)
            {
                parameters = parameters_;
                ParametersLexLSE lexlse_parameters;

                lexlse_parameters.tol_linear_dependence          = parameters.tol_linear_dependence;
                lexlse_parameters.regularization_type            = parameters.regularization_type;
                lexlse_parameters.max_number_of_CG_iterations    = parameters.max_number_of_CG_iterations;
                lexlse_parameters.variable_regularization_factor = parameters.variable_regularization_factor;

                lexlse.setParameters(lexlse_parameters);

                if (parameters.cycling_handling_enabled)
                {
                    cycling_handler.set_max_counter(parameters.cycling_max_counter);
                    cycling_handler.set_relax_step(parameters.cycling_relax_step);
                }
            }


            /**
                \brief Set data of objective ObjIndex (ObjType = GENERAL_OBJECTIVE is assumed)

                \param[in] ObjIndex Index of objective
                \param[in] data     [A,LowerBounds,UpperBounds]
            */
            void setData(Index ObjIndex, const dMatrixType& data)
            {
                if (ObjIndex >= nObj)
                {
                    throw Exception("ObjIndex >= nObj");
                }

                if (objectives[ObjIndex].getObjType() != GENERAL_OBJECTIVE)
                {
                    throw Exception("ObjType = GENERAL_OBJECTIVE is assumed");
                }

                if (objectives[ObjIndex].getDim() != data.rows())
                {
                    throw Exception("Incorrect number of equations");
                }

                // check bounds
                RealScalar bl, bu;
                for (Index CtrIndex=0; CtrIndex<objectives[ObjIndex].getDim(); CtrIndex++)
                {
                    bl = data.coeffRef(CtrIndex,nVar);
                    bu = data.coeffRef(CtrIndex,nVar+1);

                    if (isEqual(bl,bu))
                    {
                        // don't activate meaningless constraints
                        if (data.row(CtrIndex).head(nVar).squaredNorm() > 0)
                        {
                            activate(ObjIndex,CtrIndex,CTR_ACTIVE_EQ,false);
                        }
                        else
                        {
                            //printf("WARNING: equality constraint (0*x = b) not activated (obj_index = %d, ctr_index = %d) \n", ObjIndex, CtrIndex);
                        }
                    }
                    else if (bl > bu)
                    {
                        /**
                         * \todo [AS] Minor inconsistency in the API.
                         *
                         * \verbatim
                         *
                         * This results in a failure
                         *    1 <= x <= -1
                         *
                         * This is ok
                         *    1       <=  x  <=  10^10
                         *    -10^10  <=  x  <=  -1
                         *
                         * Even though the problems are equivalent.
                         *
                         * Dimitar:
                         *  We explicitly assume that upper bounds are greater
                         *  or equal to lower bounds.
                         *
                         *  It is interesting to mention as well the reason for
                         *  this assumption. Note that we use a primal
                         *  active-set method and it requires a feasible
                         *  initial iterate. There is no scalar “v” for which
                         *
                         *  1 <= x - v <= -1
                         *
                         *  is satisfied. You would need to split things in two
                         *  inequalities (in the way that you do) so that you
                         *  can play with two scalars “v_1” and “v_2”. Of
                         *  course, one could detect cases with lb > ub and
                         *  reformulate the problem but I don’t see a good
                         *  reason to do that, as I think that lb > ub simply
                         *  indicates an error in the problem formulation (do
                         *  you have an example of a meaningful problem for
                         *  which lb > ub?).
                         *
                         *  Maybe the error
                         *  "Lower bound is greater than upper bound"
                         *  should be changed so that the user clearly
                         *  understands what we assume. While the above note
                         *  could be included somewhere else.
                         * \endverbatim
                         */
                        throw Exception("(general) Lower bound is greater than upper bound.");
                    }
                }

                objectives[ObjIndex].setData(data);
            }

            /**
                \brief Set data of objective ObjIndex (ObjType = SIMPLE_BOUNDS_OBJECTIVE is assumed)

                \param[in] ObjIndex Index of objective
                \param[in] VarIndex Index variables subject to simple bounds
                \param[in] data     [LowerBounds,UpperBounds]
            */
            void setData(Index ObjIndex, Index *VarIndex, const dMatrixType& data)
            {
                if (ObjIndex >= nObj)
                {
                    throw Exception("ObjIndex >= nObj");
                }

                if (objectives[ObjIndex].getObjType() != SIMPLE_BOUNDS_OBJECTIVE)
                {
                    throw Exception("ObjType = SIMPLE_BOUNDS_OBJECTIVE is assumed");
                }

                if (objectives[ObjIndex].getDim() != data.rows())
                {
                    throw Exception("Incorrect number of equations");
                }

                // check bounds
                RealScalar bl, bu;
                for (Index CtrIndex=0; CtrIndex<objectives[ObjIndex].getDim(); CtrIndex++)
                {
                    bl = data.coeffRef(CtrIndex,0);
                    bu = data.coeffRef(CtrIndex,1);

                    if (isEqual(bl,bu))
                    {
                        activate(ObjIndex,CtrIndex,CTR_ACTIVE_EQ,false);
                    }
                    else if (bl > bu)
                    {
                        throw Exception("(simple) Lower bound is greater than upper bound.");
                    }
                }

                // check whether VarIndex contains repeated indexes (VarIndex is not assumed to be sorted)
                for (Index k=0; k<objectives[ObjIndex].getDim(); k++)
                {
                    for (Index j=0; j<objectives[ObjIndex].getDim(); j++)
                    {
                        if ((VarIndex[k] == VarIndex[j]) && (j != k))
                        {
                            throw Exception("Elements of VarIndex are not unique.");
                        }
                    }
                }

                objectives[ObjIndex].setData(VarIndex, data);
            }

            /**
                \brief Set (a non-negative) regularization factor for objective ObjIndex

                \note Regularization of an objective of type SIMPLE_BOUNDS_OBJECTIVE is not performed
            */
            void setRegularizationFactor(Index ObjIndex, RealScalar factor)
            {
                // @todo: check whether ObjIndex and factor make sense.

                objectives[ObjIndex].setRegularization(factor);
            }

            /**
                \brief Return the (primal) solution vector
            */
            dVectorType& get_x()
            {
                return x;
            }

            /**
                \brief Return the solution vector of the lexlse problem corresponding to the last
                active set

                \todo In getLambda as well I solve a lexlse problem (this is wasteful).
            */
            dVectorType& get_xStar()
            {
                formLexLSE();
                lexlse.factorize();
                lexlse.solve();

                return lexlse.get_x();
            }

            dVectorType& get_v(Index ObjIndex)
            {
                return objectives[ObjIndex].get_v();
            }

            /**
               \brief Returns the (minimal) constraint violation for objective ObjIndex

               \param[in]  ObjIndex      objective index
               \param[out] ctr_violation vector of constraint violations

               \note The result might be different from get_v() if the active-set iterations are
               prematurely terminated.
            */
            void getConstraintViolation(Index ObjIndex, dVectorType &ctr_violation)
            {
                objectives[ObjIndex].getConstraintViolation(ctr_violation);
            }

            /**
                \brief Outputs the Lagrange multipliers associated to all constraintes involved in all objectives

                \note The order of constraints is like the one provided by the user (in the problem definition)
            */
            void getLambda(std::vector<dMatrixType> & vec_lambda)
            {
                Index nActiveCtr = 0; // number of active constraints
                vec_lambda.resize(nObj);
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++) // Objectives of LexLSI
                {
                    nActiveCtr += objectives[ObjIndex].getActiveCtrCount();
                    vec_lambda[ObjIndex].setZero(getObjDim(ObjIndex),nObj);
                }

                // make sure that objectives[ObjIndex].getActiveCtrCount() is the same as
                // lexlse.getDim(ObjIndex) in case the problem is not solved.
                if (status != PROBLEM_SOLVED)
                {
                    //printf("Warning: status = %d, solving lexlse problem in getLambda(...) \n", status);
                    formLexLSE();
                    lexlse.factorize();
                }

                // "L_active" contains only the Lagrange multipliers associated to the active
                // constraints in the working set (the order in the working set is preserved).
                dMatrixType L_active = dMatrixType::Zero(nActiveCtr,nObj);
                Index nMeaningful = lexlse.getFixedVariablesCount();

                Index CtrIndex2Remove;
                int ObjIndex2Remove;
                RealScalar maxAbsValue;
                for (Index ObjIndex=0; ObjIndex<nObj-nObjOffset; ObjIndex++) // Objectives of LexLSE
                {

                    //lexlse.ObjectiveSensitivity(ObjIndex);

                    // test with the function that is actually used within the ative-set method
                    lexlse.ObjectiveSensitivity(ObjIndex,
                                                CtrIndex2Remove, ObjIndex2Remove,
                                                parameters.tol_wrong_sign_lambda,
                                                parameters.tol_correct_sign_lambda,
                                                maxAbsValue);

                    nMeaningful += lexlse.getDim(ObjIndex);
                    L_active.col(nObjOffset + ObjIndex).head(nMeaningful) = lexlse.getWorkspace().head(nMeaningful);
                }

                Index ind;
                Index accumulate_active_ctr = 0;

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++) // Objectives of LexLSI
                {
                    for (Index k=0; k<objectives[ObjIndex].getActiveCtrCount(); k++)
                    {
                        ind = objectives[ObjIndex].getActiveCtrIndex(k);
                        vec_lambda[ObjIndex].row(ind) = L_active.row(accumulate_active_ctr+k);
                    }
                    accumulate_active_ctr += objectives[ObjIndex].getActiveCtrCount();
                }
            }

            dMatrixType get_lexqr()
            {
                return lexlse.get_lexqr();
            }

            dMatrixType get_data()
            {
                return lexlse.get_data();
            }

            dMatrixType get_X_mu()
            {
                return lexlse.get_X_mu();
            }

            dMatrixType get_X_mu_rhs()
            {
                return lexlse.get_X_mu_rhs();
            }

            dVectorType get_residual_mu()
            {
                return lexlse.get_residual_mu();
            }


            /**
                \brief Get number of cycling relaxations
            */
            Index getCyclingCounter() const
            {
                return cycling_handler.get_counter();
            }

            /**
                \brief Returns number of iterations in the active-set method
            */
            Index getFactorizationsCount() const
            {
                return nFactorizations;
            }

            /**
                \brief Returns number of iterations during which a constraint has been added to the
                working set
            */
            Index getActivationsCount() const
            {
                return nActivations;
            }

            /**
                \brief Returns number of iterations during which a constraint has been removed to the
                working set
            */
            Index getDeactivationsCount() const
            {
                return nDeactivations;
            }

            Index getActiveCtrCount(Index ObjIndex) const
            {
                return objectives[ObjIndex].getActiveCtrCount();
            }

            /**
                \brief Returns number of active constraints
            */
            Index getActiveCtrCount() const
            {
                Index n = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    n += objectives[ObjIndex].getActiveCtrCount();
                }

                return n;
            }

            /**
                \brief Outputs the types (CTR_INACTIVE, CTR_ACTIVE_LB, CTR_ACTIVE_UB) of constraints for a given objective
            */
            void getActiveCtr(Index ObjIndex, std::vector<ConstraintActivationType>& ctr_type) const
            {
                Index ind;
                Index dim = objectives[ObjIndex].getDim();
                ctr_type.resize(dim,CTR_INACTIVE);
                for (Index k=0; k<objectives[ObjIndex].getActiveCtrCount(); k++)
                {
                    ind = objectives[ObjIndex].getActiveCtrIndex(k);
                    ctr_type[ind] = objectives[ObjIndex].getActiveCtrType(k);
                }
            }

            /**
                \brief Outputs the indexes and types of active constraints
            */
            void getActiveCtr_order(std::vector<ConstraintIdentifier>& ctr) const
            {
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    for (Index k=0; k<objectives[ObjIndex].getActiveCtrCount(); k++)
                    {
                        ConstraintIdentifier aCtr;
                        aCtr.set(ObjIndex,
                                 objectives[ObjIndex].getActiveCtrIndex(k),
                                 objectives[ObjIndex].getActiveCtrType(k));

                        ctr.push_back(aCtr);
                    }
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
                return objectives[ObjIndex].getDim();
            }

            /**
                \brief Returns working_set_log
            */
            std::vector<WorkingSetLogEntry>& getWorkingSetLog()
            {
                return working_set_log;
            }

            /**
               \brief Return the TotalRank of a lexlse problem
            */
            Index getTotalRank()
            {
                return lexlse.getTotalRank();
            }

        private:

            /**
                \brief Some tests on the validity of hot-start (currently only related to advanced initialization)

                \todo Additional tests shouls be implemented (e.g., feasibility of (x0,v0)).
            */
            void hot_start_related_tests()
            {
                // make sure that v0 is not only partially specified
                bool v0_is_only_partially_specified = false;
                bool user_attempted_to_spevify_v0 = objectives[0].getFlag_v0_is_specified();
                for (Index ObjIndex=1; ObjIndex<nObj; ObjIndex++)
                {
                    if (objectives[ObjIndex].getFlag_v0_is_specified() != user_attempted_to_spevify_v0)
                    {
                        // here just print a warning (actually disregard the user specified v0 below)
                        printf("WARNING: disregarding v0 because it is only partially initialized. \n");

                        user_attempted_to_spevify_v0 = true;
                        v0_is_only_partially_specified = true;
                        break;
                    }
                }

                // make sure that the user did not attempt to hot-start with {W0,v0} or {v0}
                bool user_attempted_to_spevify_v0_but_forgot_x_guess = false;
                if ((!x_guess_is_specified) && user_attempted_to_spevify_v0)
                {
                    // here just print a warning (actually disregard the user specified v0 below)
                    printf("WARNING: disregarding v0 because x_guess is not set. \n");
                    user_attempted_to_spevify_v0_but_forgot_x_guess = true;
                }

                if (v0_is_only_partially_specified || user_attempted_to_spevify_v0_but_forgot_x_guess)
                {
                    // disregard user input for v0
                    for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    {
                        objectives[ObjIndex].setFlag_v0_is_specified(false);
                    }
                }
            }

            /**
               \brief Computes an initial feasible pair (x,w)

                \verbatim
                ----------------------------------------------------------------------------------------------------
                two cases:
                ----------------------------------------------------------------------------------------------------
                1. x_guess is not specified

                    x = x_star        ,  v{active} = A{active}*x-b{active},  v{inactive} = middle of bounds
                   dx = x_star - x = 0, dv{active} = A{active}*dx = 0     , dv{inactive} = -v{inactive}

                2. x_guess is specified

                    x = x_guess       ,  v{active} = A{active}*x-b{active},  v{inactive} = middle of bounds
                   dx = x_star - x    , dv{active} = A{active}*dx         , dv{inactive} = -v{inactive}

                   dv{active} = A{active}*x_star - b{active} - (A{active}*x - b{active}) = A{active}*(x_star - x)
                ----------------------------------------------------------------------------------------------------
                \endverbatim
            */
            void phase1()
            {
                hot_start_related_tests();

                if (!x_guess_is_specified)
                {
                    formLexLSE();
                    lexlse.factorize();
                    lexlse.solve(); // solve lexlse based on initial working set
                    lexlse_rank = getTotalRank();

                    x = lexlse.get_x();
                }

                // --------------------------------------------------------
                // form initial working set and a feasible pair (x,v)
                // --------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].phase1(x, // either x = x_guess or x = lexlse.get_x()
                                                x_guess_is_specified,
                                                parameters.modify_type_active_enabled,
                                                parameters.modify_type_inactive_enabled,
                                                parameters.modify_x_guess_enabled,
                                                parameters.set_min_init_ctr_violation,
                                                parameters.tol_feasibility);
                }

                // --------------------------------------------------------
                // form dx
                // --------------------------------------------------------
                if (x_guess_is_specified)
                {
                    formLexLSE();
                    lexlse.factorize();
                    lexlse.solve(); // solve lexlse based on initial working set
                    lexlse_rank = getTotalRank();

                    dx = lexlse.get_x() - x; // take a step towards x_star
                }
                else
                {
                    // dx.setZero(); // already set to zero in initialize()
                }

                // --------------------------------------------------------
                // form step for v (similar to formStep() but dx is initialized above)
                // --------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].formStep(dx);
                }
                // --------------------------------------------------------

                nFactorizations++; // one factorization is performed
            }

            /**
               \brief An alternative to #phase1()

               \note dx = 0, dv{active} = 0, dv{inactive} = -v{inactive}

               \note A factorization is NOT performed here

               \attention x_guess has to be specified by the user.
            */
            void phase1_v0()
            {
                if (!x_guess_is_specified)
                {
                    throw Exception("when use_phase1_v0 = true, x_guess has to be specified");
                }

                hot_start_related_tests();

                // --------------------------------------------------------
                // form initial working set and a feasible pair (x,v)
                // --------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].phase1(x,                    // x = x_guess
                                                x_guess_is_specified, // true
                                                parameters.modify_type_active_enabled,
                                                parameters.modify_type_inactive_enabled,
                                                parameters.modify_x_guess_enabled,
                                                parameters.set_min_init_ctr_violation,
                                                parameters.tol_feasibility);
                }

                // --------------------------------------------------------
                // form dx
                // --------------------------------------------------------
                // dx.setZero(); // already set to zero in initialize()

                // --------------------------------------------------------
                // form step for v (similar to formStep() but dx is initialized above)
                // --------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].formStep(dx);
                }
                // --------------------------------------------------------
            }

            /**
                \brief Resize LexLSI problem

                \param[in] ObjDim_  Number of constraints involved in each objective
                \param[in] ObjType_ Type of each objective
            */
            void resize(Index *ObjDim_, ObjectiveType *ObjType_)
            {
                nObjOffset = 0;
                if (ObjType_[0] == SIMPLE_BOUNDS_OBJECTIVE) // If only simple bounds in first objective of LexLSI
                {
                    nObjOffset = 1;
                }

                // In LexLSE, fixed variables are handled separately and are not defined as an objective
                // ObjDim_ + nObjOffset is pointer arithmetic
                lexlse.resize(nVar, nObj - nObjOffset, ObjDim_ + nObjOffset);

                nActive.resize(nObj);
                objectives.resize(nObj);
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].resize(ObjDim_[ObjIndex],nVar,ObjType_[ObjIndex]);
                }

                x.resize(nVar);
                dx.resize(nVar);

                initialize();
            }

            /**
                \brief Initializations
            */
            void initialize()
            {
                nIterations     = 0;
                nActivations    = 0;
                nDeactivations  = 0;
                nFactorizations = 0;
                lexlse_rank     = 0;

                step_length = 0;

                x.setZero();
                dx.setZero();
            }

            /**
                \brief Form an LexLSE problem (using the current working set)
            */
            void formLexLSE()
            {
                // obj_info.FirstRowIndex has to be initialized before I start setting CtrType in formLexLSE below
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nActive(ObjIndex) = objectives[ObjIndex].getActiveCtrCount();
                }
                lexlse.setObjDim(&nActive(0)+nObjOffset);

                Index counter = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].formLexLSE(lexlse, counter, ObjIndex-nObjOffset);
                }
            }

            /**
               \brief Form the step (dx,dw) from the current iterate and compute the step length StepLength
            */
            void formStep()
            {
                dx = lexlse.get_x() - x;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    objectives[ObjIndex].formStep(dx);
                }
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
                                          ConstraintActivationType &CtrTypeBlocking,
                                          RealScalar &alpha)
            {
                alpha = 1;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    if (objectives[ObjIndex].checkBlockingConstraints(CtrIndexBlocking, CtrTypeBlocking, alpha, parameters.tol_feasibility))
                    {
                        ObjIndexBlocking = ObjIndex;
                    }
                }

                if (alpha < 1)
                {
                    return true; // there are blocking constraints
                }
                else
                {
                    return false;
                }
            }

            /**
               \note Probably this could be done in O(n*log(n))
            */
            Index findFirstCtrWrongSign(std::vector<ConstraintInfo> &ctr_wrong_sign)
            {
                std::vector<ConstraintInfo>::iterator it = ctr_wrong_sign.end();

                Index k = 0;
                while (it == ctr_wrong_sign.end())
                {
                    it = std::find(ctr_wrong_sign.begin(), ctr_wrong_sign.end(), WS[k]);
                    k++;
                }

                return --k;
            }

            bool findActiveCtr2Remove(Index &ObjIndex2Remove, Index &CtrIndex2Remove, RealScalar &lambda_wrong_sign)
            {
                if (parameters.deactivate_first_wrong_sign)
                {
                    return findActiveCtr2Remove_first(ObjIndex2Remove, CtrIndex2Remove, lambda_wrong_sign);
                }
                else
                {
                    return findActiveCtr2Remove_largest(ObjIndex2Remove, CtrIndex2Remove, lambda_wrong_sign);
                }
            }

            // remove first ctr with Lambda with wrong sign
            bool findActiveCtr2Remove_first(Index &ObjIndex2Remove,
                                            Index &CtrIndex2Remove,
                                            RealScalar &lambda_wrong_sign)
            {
                std::vector<ConstraintInfo> ctr_wrong_sign;

                lambda_wrong_sign = 0; // this doesn't matter (todo modify later)

                bool DescentDirectionExists = false;
                for (Index ObjIndex=0; ObjIndex<nObj-nObjOffset; ObjIndex++) // loop over objectives of LexLSE problem
                {
                    lexlse.ObjectiveSensitivity(ObjIndex,
                                                parameters.tol_wrong_sign_lambda,
                                                parameters.tol_correct_sign_lambda,
                                                ctr_wrong_sign);

                    if (ctr_wrong_sign.size() > 0)
                    {
                        DescentDirectionExists = true;
                        break;
                    }
                }

                if (DescentDirectionExists)
                {
                    for (Index k=0; k<ctr_wrong_sign.size(); k++)
                    {
                        ctr_wrong_sign[k].increment_obj_index(nObjOffset);

                        int obj_tmp = ctr_wrong_sign[k].get_obj_index();
                        int ctr_tmp = ctr_wrong_sign[k].get_ctr_index();

                        ctr_tmp = objectives[obj_tmp].getActiveCtrIndex(ctr_tmp);
                        ctr_wrong_sign[k].set_ctr_index(ctr_tmp);
                    }

                    Index k = findFirstCtrWrongSign(ctr_wrong_sign);
                    ObjIndex2Remove = (Index)WS[k].get_obj_index();
                    CtrIndex2Remove = (Index)WS[k].get_ctr_index();
                    CtrIndex2Remove = objectives[ObjIndex2Remove].getCtrIndex(CtrIndex2Remove);
                }

                return DescentDirectionExists;
            }


            /**
               \brief Finds active constraints that should be removed from the working set

               \param[out] ObjIndex2Remove Index of objective from which to remove a constraint
               \param[out] CtrIndex2Remove Index of constraint in the working set of objective ObjIndex2Remove to remove

               \return true if there are constraints to remove
            */
            bool findActiveCtr2Remove_largest(Index &ObjIndex2Remove, Index &CtrIndex2Remove, RealScalar &lambda_wrong_sign)
            {
                bool DescentDirectionExists = false;
                int ObjIndex2Remove_int;
                for (Index ObjIndex=0; ObjIndex<nObj-nObjOffset; ObjIndex++) // loop over objectives of LexLSE problem
                {
                    DescentDirectionExists = lexlse.ObjectiveSensitivity(ObjIndex,
                                                                         CtrIndex2Remove,
                                                                         ObjIndex2Remove_int,
                                                                         parameters.tol_wrong_sign_lambda,
                                                                         parameters.tol_correct_sign_lambda,
                                                                         lambda_wrong_sign);

                    if (DescentDirectionExists)
                    {
                        break;
                    }
                }

                // Note that when the first objective of LexLSI is of type SIMPLE_BOUNDS_OBJECTIVE,
                // and if a constraint is to be removed from it, ObjIndex2Remove_int = -1 (see end of LexLSE.ObjectiveSensitivity(...)).
                ObjIndex2Remove = ObjIndex2Remove_int + nObjOffset; // objective of LexLSI problem

                return DescentDirectionExists;
            }

            /**
               \brief One iteration of an active-set method
            */
            OperationType verifyWorkingSet()
            {
                // ----------------------------------------------------------------------
                Index ObjIndex2Manipulate = 0, CtrIndex2Manipulate = 0; // initialize so that the compiler doesn't complain
                ConstraintActivationType CtrType2Manipulate = CTR_INACTIVE;

                bool normalIteration = true;
                OperationType operation = OPERATION_UNDEFINED;
                ConstraintIdentifier constraint_identifier(0,0,CTR_INACTIVE,0); // initialize so that g++ does not complain

                RealScalar alpha;

                bool cycling_detected;
                // ----------------------------------------------------------------------

                if (nIterations != 0) // nIterations == 0 is handled in phase1()
                {
                    formLexLSE();

                    lexlse.factorize();
                    lexlse.solve();
                    lexlse_rank = getTotalRank();

                    formStep();

                    nFactorizations++;
                }
                else // if this is the first time we enter verifyWorkingSet()
                {
                    if (parameters.use_phase1_v0) // if we have used phase1_v0()
                    {
                        normalIteration = false; // i.e., only check for blocking constraints and make a step
                    }
                }

                if (checkBlockingConstraints(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate, alpha))
                {
                    if (parameters.cycling_handling_enabled)
                    {
                        constraint_identifier.set(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate);
                    }

                    if (parameters.log_working_set_enabled)
                    {
                        WorkingSetLogEntry wlog(ObjIndex2Manipulate,
                                                CtrIndex2Manipulate,
                                                CtrType2Manipulate,
                                                alpha,
                                                lexlse_rank);

                        working_set_log.push_back(wlog);
                    }

                    operation = OPERATION_ADD;
                    activate(ObjIndex2Manipulate, CtrIndex2Manipulate, CtrType2Manipulate);
                }
                else
                {
                    if (normalIteration)
                    {
                        RealScalar lambda_wrong_sign;
                        if (findActiveCtr2Remove(ObjIndex2Manipulate, CtrIndex2Manipulate, lambda_wrong_sign))
                        {
                            if (parameters.cycling_handling_enabled)
                            {
                                constraint_identifier.set(ObjIndex2Manipulate,
                                                          objectives[ObjIndex2Manipulate].getActiveCtrIndex(CtrIndex2Manipulate),
                                                          objectives[ObjIndex2Manipulate].getActiveCtrType(CtrIndex2Manipulate));
                            }

                            if (parameters.log_working_set_enabled)
                            {
                                WorkingSetLogEntry wlog(ObjIndex2Manipulate,
                                                        objectives[ObjIndex2Manipulate].getActiveCtrIndex(CtrIndex2Manipulate),
                                                        CTR_INACTIVE,
                                                        lambda_wrong_sign,
                                                        lexlse_rank);

                                working_set_log.push_back(wlog);
                            }

                            operation = OPERATION_REMOVE;
                            deactivate(ObjIndex2Manipulate, CtrIndex2Manipulate);
                        }
                        else
                        {
                            status = PROBLEM_SOLVED;
                        }
                    }
                }

                if (operation == OPERATION_ADD)
                {
                    step_length = alpha; // record the value of alpha
                }
                else
                {
                    step_length = -1; // this is used only for debugging purposes
                }

                if (alpha > 0) // take a step
                {
                    x += alpha*dx;
                    for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                    {
                        objectives[ObjIndex].step(alpha);
                    }
                }

                if (parameters.cycling_handling_enabled && operation != OPERATION_UNDEFINED)
                {
                    status = cycling_handler.update(operation,
                                                    constraint_identifier,
                                                    objectives,
                                                    cycling_detected);

                    if (parameters.log_working_set_enabled)
                    {
                        working_set_log.back().cycling_detected = cycling_detected;
                    }
                }

                nIterations++;

                return operation;
            }

            /**
               \brief Output stuff

               \note this file makes sence only when using phase1() (and not phase1_v0())
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

                file << "% ==============================================\n";
                file << "% nIterations       = " << nIterations << "\n";
                file << "% status            = " << status << "\n";
                file << "% counter (cycling) = " << getCyclingCounter() << "\n";
                file << "nFactorizations_("<<nIterations+1<<") = " << getFactorizationsCount() << ";\n";
                if (nIterations != 0)
                {
                    file << "operation_("<<nIterations+1<<")       = " << operation << ";\n";
                    file << "stepLength_("<<nIterations+1<<")      = " << step_length << ";\n";
                }
                file << "% ==============================================\n";

                dVectorType xStar = lexlse.get_x();

                file << "% obtained with old working set" << "\n";
                file << "xStar_(:,"<<nIterations+1<<") = [ ";
                for (Index k=0; k<nVar; k++)
                {
                    file << xStar(k) << " ";
                }
                file << "]'; \n";

                file << "% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

                file << "dx_(:,"<<nIterations+1<<") = [ ";
                for (Index k=0; k<nVar; k++)
                {
                    file << dx(k) << " ";
                }
                file << "]'; \n";

                file << "% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    dVectorType dw_ = objectives[ObjIndex].get_dv();

                    file << "dw_{"<<ObjIndex+1<<"}(:,"<<nIterations+1<<") = [ ";
                    for (Index k=0; k<objectives[ObjIndex].getDim(); k++)
                    {
                        file << dw_(k) << " ";
                    }
                    file << "]';\n";
                }

                file << "% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

                file << "x_(:,"<<nIterations+1<<") = [ ";
                for (Index k=0; k<nVar; k++)
                    file << x(k) << " ";
                file << "]'; \n";

                file << "% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    dVectorType w_ = objectives[ObjIndex].get_v();

                    file << "w_{"<<ObjIndex+1<<"}(:,"<<nIterations+1<<") = [ ";
                    for (Index k=0; k<objectives[ObjIndex].getDim(); k++)
                    {
                        file << w_(k) << " ";
                    }
                    file << "]';\n";
                }

                file << "% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    file << "a_{"<<ObjIndex+1<<"}(:,"<<nIterations+1<<") = [ ";
                    for (Index k=0; k<objectives[ObjIndex].getDim(); k++)
                    {
                        file << (Index) objectives[ObjIndex].getCtrType(k) << " ";
                    }
                    file << "]';\n";
                }

                /*
                  dMatrixType L = getLambda();

                  file << "L_{"<<nIterations+1<<"} = [";
                  for (Index i=0; i<L.rows(); i++)
                  {
                  for (Index j=0; j<L.cols(); j++)
                  {
                  file << L(i,j) << " ";
                  }
                  file << "\n";
                  }
                  file << "]; \n";
                */

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
                its ObjType = SIMPLE_BOUNDS_OBJECTIVE), the number of objectives in LexLSI and LexLSE differ
                with 1 because fixed variables are not treated as an objective in LexLSE.
            */
            Index nObjOffset;

            /*
              \brief Number of iterations during which a constraint was added
            */
            Index nActivations;

            /*
              \brief Number of iterations during which a constraint was removed
            */
            Index nDeactivations;

            /*
              \brief Number of factorization
            */
            Index nFactorizations;

            /*
              \brief Iterations counter
            */
            Index nIterations;

            /*
              \brief Stores the rank of the constraints in the last solved lexlse problem
            */
            Index lexlse_rank;

            /**
                \brief If x_guess_is_specified == true, the function set_x0(dVectorType &x0) has been
                called and x0 has been initialized.

                \note This is later used in phase1().
            */
            bool x_guess_is_specified;

            /*
              \brief Equal to alpha in verifyWorkingSet()

              \note For output/debugging purposes
            */
            RealScalar step_length;

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
            std::vector<Objective> objectives;

            /**
                \brief Vector containing activation/deactivation info for each iteration
            */
            std::vector<WorkingSetLogEntry> working_set_log;

            /**
                \brief Handle cycling
            */
            CyclingHandler cycling_handler;

            /**
                \brief Parameters of the solver.
            */
            ParametersLexLSI parameters;

            /**
                \brief Store Working Set but in order to addition (objectives are not separated)
            */
            std::vector<ConstraintInfo> WS;

        }; // END class LexLSI

    } // END namespace internal

} // END namespace LexLS

#endif // LEXLSE
