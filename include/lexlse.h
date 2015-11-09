#ifndef LEXLSE
#define LEXLSE

#include <utility.h>

namespace LexLS
{
    namespace internal
    {
        /**
            \brief Definition of a lexicographic least-squares problem with equality constraints

            \todo To describe how variables can be fixed.

            \todo Should I resize fixed_var_type.resize(nVar) and fixed_var_index.resize(nVar) only once (in
            active-set iterations)?

            \todo Now I initialize memory for many additional matrices related to the
            regularization. This could be improved.

            \todo The way we detect singularity now is not good: (maxColNormValue <
            parameters.tol_linear_dependence). We should use at least the way it is done in
            Eigen. It is easy to directly take it, but I had some problem with this in the past. We
            should have many test problem on which to evaluate various singularity detection
            approaches (and only then to start changing things like this).

            \todo to have a custom implementation of applyOnTheLeft(householderSequence (no need to
            realocate workspace every time)
        */
        class LexLSE
        {
        public:

            // =================================================================================================
            // Constructors
            // =================================================================================================

            /**
                \brief Default constructor
            */
            LexLSE():
                nVarFixed(0),
                nVarFixedInit(0){}

            /**
                \param[in] nVar_   Number of variables (only number of elements in x, and not in the residuals v)
                \param[in] nObj_   Number of objectives
                \param[in] ObjDim_ Number of constraints involved in each objective
            */
            LexLSE(Index nVar_, Index nObj_, Index *ObjDim_):
                nVarFixed(0),
                nVarFixedInit(0)
            {
                resize(nVar_, nObj_, ObjDim_);
                setObjDim(ObjDim_);
            }

            // =================================================================================================
            //
            // =================================================================================================

            /**
                \brief Allocate memory

                \param[in] nVar_     Number of variables
                \param[in] nObj_     Number of objectives
                \param[in] maxObjDim Maximum number of constraints involved in each objective
            */
            void resize(Index nVar_, Index nObj_, Index *maxObjDim)
            {
                nVar = nVar_;
                nObj = nObj_;

                obj_info.resize(nObj);
                x.resize(nVar);
                P.resize(nVar);

                Index maxObjDimSum = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    maxObjDimSum += maxObjDim[ObjIndex];
                }

                hh_scalars.resize(maxObjDimSum);
                ctr_type.resize(maxObjDimSum, CTR_INACTIVE); // used only for initialization purposes

                LOD.resize(maxObjDimSum,nVar+1); // store the RHS as well, thus the "+1"
                PROBLEM_DATA.resize(maxObjDimSum,nVar+1);

                Index dim = std::max(maxObjDimSum,nVar);
                dWorkspace.resize(2*dim + nVar + 1);

                column_permutations.resize(nVar);

                null_space.resize(nVar,nVar+1);
                array.resize(nVar,nVar+1);

                X_mu.resize(nVar,nObj);
                X_mu_rhs.resize(nVar,nObj);

                residual_mu.resize(maxObjDimSum); // initialized in factorize(...)

                // no need to initialize them (used in cg_tikhonov(...))
                rqsp_work.resize(2*nVar,4);
            }

            /**
                \brief Factorize using Householder QR (with column pivoting) + classical implementation of Gauss elimination

                \note This is an unblocked implementation of both the QR decomposition (essentially,
                Eigen's ColPivHouseholderQR) and the Gauss elimination.

                \attention When the rank of the constraints involved in a given objective is equal to 1,
                many simplifications can be made. In the current implementation, this case is benefited
                from only when forming the Schur complement (during the Gauss elimination).

                \todo residual_mu and fixed variables (see WARNING below)
            */
            void factorize()
            {
                PROBLEM_DATA = LOD;

                // ================================================================================
                // ================================================================================
                // ================================================================================
                if (false) // todo: remove after testing
                {
                    aRegularizationFactor = 0.0;

                    std::cout << "nVar                  = " << nVar                  << "\n"
                              << "nObj                  = " << nObj                  << "\n"
                              << "nCtr                  = " << nCtr                  << "\n"
                              << "nVarFixed             = " << nVarFixed             << "\n"
                              << "nVarFixedInit         = " << nVarFixedInit         << "\n"
                              << "TotalRank             = " << TotalRank             << "\n"
                              << "aRegularizationFactor = " << aRegularizationFactor << "\n"
                              << "fixed_var_type.size() = " << fixed_var_type.size() << "\n";

                    for (Index kk=0; kk<nObj; kk++)
                        obj_info[kk].print();
                    parameters.print();

                    std::cout << "[";
                    for (Index kk=0; kk<nCtr; kk++)
                        std::cout << ctr_type[kk] << " ";
                    std::cout << "]\n";

                    fixed_var_index.setZero();
                    P.setIdentity();
                    column_permutations.setZero();
                    dWorkspace.setZero();
                    x.setZero();
                    hh_scalars.setZero();
                    residual_mu.setZero();
                    null_space.setZero();
                    array.setZero();
                    rqsp_work.setZero();
                    X_mu.setZero();
                    X_mu_rhs.setZero();
                    std::cout << "-----------------------------------------------------\n\n";
                }
                // ================================================================================
                // ================================================================================
                // ================================================================================


                RealScalar maxColNormValue, tau, PivotValue;
                Index RemainingRows, ObjRank, ObjDim, maxColNormIndex;

                Index RowIndex;                   // Current constraint
                Index FirstRowIndex;              // The same as obj_info[k].first_row_index (already available)
                Index FirstColIndex;              // The same as obj_info[k].first_col_index (not yet computed)
                Index FirstRowIndexNextObjective; // The first index of the next objective

                // --------------------------------------------------------------------------
                // Handling of fixed valiables (apply row & column permutations)
                // --------------------------------------------------------------------------
                if (nVarFixed>0) // if there are fixed variables
                {
                    Index coeff;

                    // Permute variables so that the fixed variables come first (FIXME: explain better)
                    for (Index k=0; k<nVarFixed; k++)
                    {
                        coeff = fixed_var_index.coeffRef(k);
                        column_permutations.coeffRef(k) = coeff;
                        if (k != coeff)
                        {
                            LOD.col(k).head(nCtr).swap(LOD.col(coeff).head(nCtr));
                        }

                        for (Index i=k+1; i<nVarFixed; i++)
                        {
                            if (fixed_var_index.coeffRef(i) == k)
                            {
                                fixed_var_index.coeffRef(i) = coeff; // fixed_var_index is modified, but this is OK.
                                break;
                            }
                        }
                    }
                    LOD.col(nVar).head(nCtr).noalias() -= LOD.block(0,0,nCtr,nVarFixed)*x.head(nVarFixed);
                }
                // --------------------------------------------------------------------------

                // Jump over the first nVarFixed columns (they are dedicated to the fixed variables)
                // There is no jump in the rows because the identity matrix is not stored explicitly
                Index ColIndex         = nVarFixed;        // Current variable
                Index RemainingColumns = nVar - nVarFixed; // Remove the fixed variables from the available variables

                if (ColIndex >= nVar)
                {
                    TotalRank = nVarFixed;

                    // form the permutation matrix
                    for (Index k=0; k<TotalRank; k++)
                    {
                        P.applyTranspositionOnTheRight(k, column_permutations.coeff(k));
                    }

                    return; // early termination if possible
                }

                // ----------------------------------------------------
                dVectorBlockType    ColNorms(dWorkspace,    0, nVar  );
                dVectorBlockType hhWorkspace(dWorkspace, nVar, nVar+1);
                // ----------------------------------------------------

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++) // loop over all (explicitly defined) objectives
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index = ColIndex;
                    ObjDim        = obj_info[ObjIndex].dim;

                    // residual_mu.segment(FirstRowIndex,ObjDim) has to be initialized before the
                    // Householder transformations but after the elimination steps.
                    // WARNING: how about fixed variables (no time to think about this now)
                    residual_mu.segment(FirstRowIndex,ObjDim) = LOD.col(nVar).segment(FirstRowIndex,ObjDim);

                    for(Index k=ColIndex; k<nVar; k++) // initially compute the norms of the columns
                    {
                        ColNorms.coeffRef(k) = LOD.col(k).segment(FirstRowIndex,ObjDim).squaredNorm();
                    }

                    // QR factorization of constraints involved in objective ObjIndex using the remaining variables
                    for (Index counter=0; counter<ObjDim; counter++) // loop over all constraints in a given objective
                    {
                        RowIndex      = FirstRowIndex + counter; // current row to process
                        RemainingRows = ObjDim        - counter; // remaining rows in the current objective

                        maxColNormValue = ColNorms.tail(RemainingColumns).maxCoeff(&maxColNormIndex);
                        maxColNormIndex += ColIndex;

                        // the next two lines are for numerical stability
                        // (I use them because this is what is done in ColPivHouseholderQR.h)
                        maxColNormValue = LOD.col(maxColNormIndex).segment(RowIndex,RemainingRows).squaredNorm();
                        ColNorms.coeffRef(maxColNormIndex) = maxColNormValue;

                        // After we break, elimination is performed if there are objectives with lower priority and obj_info[ObjIndex].rank > 0
                        if (maxColNormValue < parameters.tol_linear_dependence)
                        {
                            break;
                        }

                        // --------------------------------------------------------------------------
                        // apply column permutation
                        // --------------------------------------------------------------------------
                        column_permutations.coeffRef(ColIndex) = maxColNormIndex;
                        if(ColIndex != maxColNormIndex)
                        {
                            LOD.col(ColIndex).head(nCtr).swap(LOD.col(maxColNormIndex).head(nCtr));
                            std::swap(ColNorms.coeffRef(ColIndex), ColNorms.coeffRef(maxColNormIndex));

                            // FOR REGULARIZATION
                            null_space.col(ColIndex)
                                .head(FirstColIndex).swap(null_space.col(maxColNormIndex).head(FirstColIndex));
                        }

                        // --------------------------------------------------------------------------
                        // apply Householder transformations (on the RHS as well)
                        // --------------------------------------------------------------------------
                        // when RemainingRows = 1, since sqrt(maxColNormValue) >= parameters.tol_linear_dependence,
                        // the Householder matrix is the identity (tau = 0)
                        if (RemainingRows > 1)
                        {
                            LOD.col(ColIndex).segment(RowIndex,RemainingRows).makeHouseholderInPlace(tau,PivotValue);
                            LOD.coeffRef(RowIndex,ColIndex) = PivotValue;
                            LOD.block(RowIndex,ColIndex+1,RemainingRows,RemainingColumns) // apply transformation on the RHS as well
                                .applyHouseholderOnTheLeft(LOD.col(ColIndex).segment(RowIndex+1,RemainingRows-1),
                                                           tau,
                                                           &hhWorkspace.coeffRef(0));
                            hh_scalars.coeffRef(FirstRowIndex+counter) = tau;
                        }
                        // --------------------------------------------------------------------------

                        ColIndex++;
                        RemainingColumns = nVar - ColIndex;

                        // terminate the QR factorization (after the elimination step below, the LOD is terminated as well)
                        if (RemainingColumns == 0)
                        {
                            break;
                        }

                        // update our table of squared norms of the columns
                        // (note that above ColIndex is incremented and RemainingColumns is decreased)
                        if (RemainingRows > 0)
                        {
                            ColNorms.tail(RemainingColumns) -= LOD.row(RowIndex).segment(ColIndex,RemainingColumns).cwiseAbs2();
                        }

                    } // END for (counter=0; counter<ObjDim; counter++)

                    // Note that here ColIndex is the index of the next available variable

                    ObjRank = obj_info[ObjIndex].rank = ColIndex - FirstColIndex; // store the rank

                    // -----------------------------------------------------------------------
                    // Regularization
                    // -----------------------------------------------------------------------
                    if (0.0 == parameters.variable_regularization_factor) // constant regularization factor
                    {
                        aRegularizationFactor = obj_info[ObjIndex].regularization_factor;
                    }
                    else // variable damping factor (use with caution)
                    {
                        aRegularizationFactor = 0.0;
                        if (ObjRank > 0)
                        {
                            // -----------------------------------------------------------------------
                            // conditioninig estimation
                            // -----------------------------------------------------------------------
                            dVectorType rhs_tmp = LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                            RealScalar conditioning_estimate = rhs_tmp.squaredNorm();
                            LOD.block(FirstRowIndex, FirstColIndex, ObjRank, ObjRank).
                                triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(rhs_tmp);

                            conditioning_estimate /= rhs_tmp.squaredNorm();
                            // -----------------------------------------------------------------------

                            /*
                              see equation (10) in
                              Stefano Chiaverini, Bruno Siciliano, "Review of the damped least-squares inverse kinematics
                              with experiments on an industrial robot manipulator, " 1994.
                            */
                            RealScalar epsilon = parameters.variable_regularization_factor;
                            if (conditioning_estimate < epsilon)
                            {
                                aRegularizationFactor = std::sqrt((1 - (conditioning_estimate*conditioning_estimate)/(epsilon*epsilon)));
                                aRegularizationFactor *= obj_info[ObjIndex].regularization_factor;
                            }
                        }
                    }

                    if (true) //(ObjRank > 0)
                    {
                        switch(parameters.regularization_type)
                        {
                            case REGULARIZATION_TIKHONOV:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    if (FirstColIndex + ObjRank <= RemainingColumns)
                                    {
                                        regularize_tikhonov_2(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                    }
                                    else
                                    {
                                        regularize_tikhonov_1(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                    }
                                }
                                accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

                                break;

                            case REGULARIZATION_TIKHONOV_CG:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_tikhonov_CG(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                    //regularize_tikhonov_CG_x0(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                break;

                            case REGULARIZATION_R:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_R(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                break;

                            case REGULARIZATION_R_NO_Z:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_R_NO_Z(FirstRowIndex, FirstColIndex, ObjRank);
                                }
                                break;

                            case REGULARIZATION_RT_NO_Z:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_RT_NO_Z(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                break;

                            case REGULARIZATION_RT_NO_Z_CG:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_RT_NO_Z_CG(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                break;

                            case REGULARIZATION_TIKHONOV_1:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    //regularize_tikhonov_1(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                    regularize_tikhonov_1_test(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns, ObjIndex);
                                }
                                accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                break;

                            case REGULARIZATION_TIKHONOV_2:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_tikhonov_2(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                break;

                            case REGULARIZATION_TEST:

                                if ( !isEqual(aRegularizationFactor,0.0) )
                                {
                                    regularize_test(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                                }
                                break;

                            case REGULARIZATION_NONE:

                                // do nothing
                                break;
                        }
                    }
                    /*
                    else
                    {
                        if (parameters.regularization_type == 7)
                        {
                            if (ObjIndex > 0)
                            {
                                X_mu.col(ObjIndex) = X_mu.col(ObjIndex-1);
                                X_mu_rhs.col(ObjIndex) = X_mu_rhs.col(ObjIndex-1);
                                residual_mu.segment(FirstRowIndex,ObjDim) = -LOD.col(nVar).segment(FirstRowIndex,ObjDim);
                            }
                        }
                    }
                    */

                    // -----------------------------------------------------------------------
                    // Gauss transformation
                    // -----------------------------------------------------------------------

                    if (ObjIndex < nObj-1) // if there are objectives with lower priority
                    {
                        if (ObjRank > 0) // if there are variables to eliminate
                        {
                            FirstRowIndexNextObjective = FirstRowIndex + ObjDim;
                            RemainingRows = nCtr-FirstRowIndexNextObjective; // VARIABLE REDEFINITION: remaining rows after the current objective

                            // update the trailing block of the matrix LOD
                            // handle the RHS vector as well, hence the "+1" (recall that RemainingColumns = nVar-ColIndex)
                            // here, we cannot use directly LOD.bottomRightCorner(RemainingRows,RemainingColumns+1).noalias()

                            dBlockType LeftBlock(LOD,FirstRowIndexNextObjective, FirstColIndex, RemainingRows, ObjRank);
                            dBlockType UpBlock(LOD,FirstRowIndex,ColIndex,ObjRank,RemainingColumns+1);
                            dBlockType TrailingBlock(LOD,FirstRowIndexNextObjective,ColIndex,RemainingRows,RemainingColumns+1);

                            LOD.block(FirstRowIndex,FirstColIndex,ObjRank,ObjRank)
                                .triangularView<Eigen::Upper>()
                                .solveInPlace<Eigen::OnTheRight>(LeftBlock);

                            if (ObjRank == 1)
                            {
                                // the .col(0) and .row(0) are important only for efficient computation
                                TrailingBlock.noalias() -= LeftBlock.col(0) * UpBlock.row(0);
                            }
                            else if (ObjRank > 2 && ObjRank <= 8)
                            {
                                for (Index k=0; k<RemainingColumns+1; k++)
                                {
                                    TrailingBlock.col(k).noalias() -= LeftBlock * UpBlock.col(k);
                                }
                            }
                            else if ((ObjRank == 2) || (ObjRank > 8))
                            {
                                TrailingBlock.noalias() -= LeftBlock * UpBlock;
                            }
                        }
                    }

                    // -----------------------------------------------------------------------

                    if (RemainingColumns == 0) // ColIndex = nVar
                    {
                        // Initialize some remaining fields of obj_info before we terminate (used later in ComputeLambda())
                        for (Index k=ObjIndex+1; k<nObj; k++)
                        {
                            // of course, obj_info[>ObjIndex].rank = 0
                            obj_info[k].first_col_index = obj_info[k-1].first_col_index + obj_info[k-1].rank;

                            X_mu.col(k) = X_mu.col(k-1);
                            X_mu_rhs.col(k) = X_mu_rhs.col(k-1);
                            residual_mu.segment(obj_info[k].first_row_index,obj_info[k].dim) = \
                                - LOD.col(nVar).segment(obj_info[k].first_row_index,obj_info[k].dim);
                        }

                        break; // terminate the LOD
                    }

                } // END for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)

                TotalRank = nVarFixed;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    TotalRank += obj_info[ObjIndex].rank;
                }

                // form the permutation matrix
                for (Index k=0; k<TotalRank; k++)
                {
                    P.applyTranspositionOnTheRight(k, column_permutations.coeff(k));
                }

            } // END factorize()

            /**
                \brief Output all constraints whose Lambda is with wrong sign
            */
            void ObjectiveSensitivity(Index ObjIndex,
                                      RealScalar tol_wrong_sign_lambda,
                                      RealScalar tol_correct_sign_lambda,
                                      std::vector<ConstraintInfo> &ctr_wrong_sign)
            {
                Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
                Index &ColDim = FirstColIndex; // ColDim will be used to indicate column dimension (I use it for clarity)

                // Even though this computation can be improved, it is very cheap and it is convenient to perform here
                Index nLambda = 0; // Number of constraints that may influence objective ObjIndex (excluding "variable fixing" constraints)
                Index nRank   = 0; // Total rank of objectives before objective ObjIndex
                for (Index k=0; k<ObjIndex; k++)
                {
                    nLambda += obj_info[k].dim;
                    nRank   += obj_info[k].rank;
                }
                nLambda += obj_info[ObjIndex].dim;

                // ---------------------------------------------------------------------------------
                dWorkspace.head(nVarFixed + nLambda + nRank + nVarFixed).setZero();
                dVectorBlockType LambdaFixed(dWorkspace,                   0, nVarFixed);
                dVectorBlockType      Lambda(dWorkspace,           nVarFixed, nLambda);
                dVectorBlockType         rhs(dWorkspace, nLambda + nVarFixed, nRank+nVarFixed);
                // ---------------------------------------------------------------------------------

                FirstRowIndex = obj_info[ObjIndex].first_row_index;
                FirstColIndex = obj_info[ObjIndex].first_col_index;
                ObjDim        = obj_info[ObjIndex].dim;
                ObjRank       = obj_info[ObjIndex].rank;

                // Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed

                // copy only what is needed to compute the residual v = A*x-b (i.e., -y_hat)
                Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = \
                    -LOD.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);

                // compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
                Lambda.segment(FirstRowIndex, ObjDim)
                    .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                  FirstColIndex,
                                                                  ObjDim,
                                                                  ObjRank),
                                                        hh_scalars.segment(FirstRowIndex,ObjDim)));

                // check for wrong sign of the Lagrange multipliers
                findDescentDirection(ObjIndex,
                                     FirstRowIndex,
                                     ObjDim,
                                     Lambda,
                                     tol_wrong_sign_lambda,
                                     tol_correct_sign_lambda,
                                     ctr_wrong_sign);

                if (ObjIndex>0) // the first objective has only Lagrange multipliers equal to the optimal residual
                {
                    // e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}
                    rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                        Lambda.segment(FirstRowIndex, ObjDim);

                    //for (int k=ObjIndex-1; k>=0; k--)
                    for (Index k=ObjIndex; k--; ) //ObjIndex-1, ..., 0.
                    {
                        FirstRowIndex = obj_info[k].first_row_index;
                        FirstColIndex = obj_info[k].first_col_index;
                        ObjDim        = obj_info[k].dim;
                        ObjRank       = obj_info[k].rank;

                        // Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed

                        Lambda.segment(FirstRowIndex, ObjRank) = rhs.segment(FirstColIndex, ObjRank);

                        // apply Q_k' on the left
                        Lambda.segment(FirstRowIndex, ObjDim)
                            .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                          FirstColIndex,
                                                                          ObjDim,
                                                                          ObjRank),
                                                                hh_scalars.segment(FirstRowIndex,ObjDim)));

                        rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                            Lambda.segment(FirstRowIndex, ObjDim);

                        // check for wrong sign of the Lagrange multiplier
                        findDescentDirection(k,
                                             FirstRowIndex,
                                             ObjDim,
                                             Lambda,
                                             tol_wrong_sign_lambda,
                                             tol_correct_sign_lambda,
                                             ctr_wrong_sign);
                    } // END for (Index k=ObjIndex-1; k>=0; k--)
                }

                if (nVarFixed>0) // Handle fixed variables (if any)
                {
                    LambdaFixed = -LOD.block(0, 0, nLambda, nVarFixed).transpose() * Lambda;

                    // check for wrong sign of the Lagrange multiplier
                    // -1(-st) objective implies the fixed variables
                    //   :if there are no fixed variables we will not be here
                    //   :if there are fixed variables ObjIndex2Remove + ObjOffset is used in LexLSI
                    findDescentDirection(-1,
                                         -1,
                                         ObjDim,
                                         Lambda,
                                         tol_wrong_sign_lambda,
                                         tol_correct_sign_lambda,
                                         ctr_wrong_sign);
                }
            }

            /**
                \brief Computes the sensitivity of objective ObjIndex with respect to (small) variatoins
                of the constraints involved in the LexLSE problem

                \note Upon exit, the lagrange multipliers associated with objective ObjIndex can be accessed using
                dWorkspace.head(nVarFixed + nLambda)
            */
            bool ObjectiveSensitivity(Index ObjIndex,
                                      Index &CtrIndex2Remove, int &ObjIndex2Remove,
                                      RealScalar tol_wrong_sign_lambda, RealScalar tol_correct_sign_lambda,
                                      RealScalar &maxAbsValue)
            {
                maxAbsValue = 0.0;
                bool tmp_bool, FoundBetterDescentDirection = false;
                Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
                Index &ColDim = FirstColIndex; // ColDim will be used to indicate column dimension (I use it for clarity)

                // Even though this computation can be improved, it is very cheap and it is convenient to perform here
                Index nLambda = 0; // Number of constraints that may influence objective ObjIndex (excluding "variable fixing" constraints)
                Index nRank   = 0; // Total rank of objectives before objective ObjIndex
                for (Index k=0; k<ObjIndex; k++)
                {
                    nLambda += obj_info[k].dim;
                    nRank   += obj_info[k].rank;
                }
                nLambda += obj_info[ObjIndex].dim;

                // ---------------------------------------------------------------------------------
                dWorkspace.head(nVarFixed + nLambda + nRank + nVarFixed).setZero();
                dVectorBlockType LambdaFixed(dWorkspace,                   0, nVarFixed);
                dVectorBlockType      Lambda(dWorkspace,           nVarFixed, nLambda);
                dVectorBlockType         rhs(dWorkspace, nLambda + nVarFixed, nRank+nVarFixed);
                // ---------------------------------------------------------------------------------

                FirstRowIndex = obj_info[ObjIndex].first_row_index;
                FirstColIndex = obj_info[ObjIndex].first_col_index;
                ObjDim        = obj_info[ObjIndex].dim;
                ObjRank       = obj_info[ObjIndex].rank;

                bool compute_residual_from_factorization = true;
                if (parameters.regularization_type == 7) // WARNING: TESTING
                {
                    compute_residual_from_factorization = false;
                }

                if (compute_residual_from_factorization)
                {
                    // Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed

                    // copy only what is needed to compute the residual v = A*x-b (i.e., -y_hat)
                    Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = \
                        -LOD.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);

                    // compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
                    Lambda.segment(FirstRowIndex, ObjDim)
                        .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                      FirstColIndex,
                                                                      ObjDim,
                                                                      ObjRank),
                                                            hh_scalars.segment(FirstRowIndex,ObjDim)));
                }
                else
                {
                    initialize_rhs(ObjIndex,rhs);

                    /*
                    // Set residual according to the constraint type. For example, if a'*x < ub
                    // (when the upper bound "ub" is active), set residual to zero. JUST TESTING.
                    RealScalar aResidual;
                    ConstraintActivationType aCtrType;
                    for (Index k=0; k<ObjDim; k++)
                    {
                        aResidual = residual_mu.segment(FirstRowIndex,ObjDim).coeff(k);
                        aCtrType  = ctr_type[obj_info[ObjIndex].first_row_index + k];
                        if ((aCtrType == CTR_ACTIVE_LB) && (aResidual > 0))
                        {
                            residual_mu.segment(FirstRowIndex,ObjDim).coeffRef(k) = 0;
                        }
                        else if ((aCtrType == CTR_ACTIVE_UB) && (aResidual < 0))
                        {
                            residual_mu.segment(FirstRowIndex,ObjDim).coeffRef(k) = 0;
                        }
                    }
                    */

                    Lambda.segment(FirstRowIndex, ObjDim) = residual_mu.segment(FirstRowIndex,ObjDim);
                }

                // check for wrong sign of the Lagrange multipliers
                FoundBetterDescentDirection = findDescentDirection(FirstRowIndex,
                                                                   ObjDim,
                                                                   maxAbsValue,
                                                                   CtrIndex2Remove,
                                                                   Lambda,
                                                                   tol_wrong_sign_lambda,
                                                                   tol_correct_sign_lambda);

                if (FoundBetterDescentDirection)
                {
                    ObjIndex2Remove = ObjIndex;
                }

                if (ObjIndex>0) // the first objective has only Lagrange multipliers equal to the optimal residual
                {
                    // e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}
                    rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                        Lambda.segment(FirstRowIndex, ObjDim);

                    //for (int k=ObjIndex-1; k>=0; k--)
                    for (Index k=ObjIndex; k--; ) //ObjIndex-1, ..., 0.
                    {
                        FirstRowIndex = obj_info[k].first_row_index;
                        FirstColIndex = obj_info[k].first_col_index;
                        ObjDim        = obj_info[k].dim;
                        ObjRank       = obj_info[k].rank;

                        // Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed

                        Lambda.segment(FirstRowIndex, ObjRank) = rhs.segment(FirstColIndex, ObjRank);

                        // apply Q_k' on the left
                        Lambda.segment(FirstRowIndex, ObjDim)
                            .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                          FirstColIndex,
                                                                          ObjDim,
                                                                          ObjRank),
                                                                hh_scalars.segment(FirstRowIndex,ObjDim)));

                        rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                            Lambda.segment(FirstRowIndex, ObjDim);

                        // check for wrong sign of the Lagrange multipliers
                        tmp_bool = findDescentDirection(FirstRowIndex,
                                                        ObjDim,
                                                        maxAbsValue,
                                                        CtrIndex2Remove,
                                                        Lambda,
                                                        tol_wrong_sign_lambda,
                                                        tol_correct_sign_lambda);

                        if (tmp_bool)
                        {
                            ObjIndex2Remove = k;
                        }
                        FoundBetterDescentDirection = (FoundBetterDescentDirection || tmp_bool);
                    } // END for (Index k=ObjIndex-1; k>=0; k--)
                }

                if (nVarFixed>0) // Handle fixed variables (if any)
                {
                    LambdaFixed = -LOD.block(0, 0, nLambda, nVarFixed).transpose() * Lambda;

                    // check for wrong sign of the Lagrange multipliers
                    tmp_bool = findDescentDirection(-1,
                                                    nVarFixed,
                                                    maxAbsValue,
                                                    CtrIndex2Remove,
                                                    LambdaFixed,
                                                    tol_wrong_sign_lambda,
                                                    tol_correct_sign_lambda);

                    // -1(-st) objective implies the fixed variables
                    //   :if there are no fixed variables we will not be here
                    //   :if there are fixed variables ObjIndex2Remove + ObjOffset is used in LexLSI
                    if (tmp_bool)
                    {
                        ObjIndex2Remove = -1;
                    }
                    FoundBetterDescentDirection = (FoundBetterDescentDirection || tmp_bool);
                }

                return FoundBetterDescentDirection;

            } // END ObjectiveSensitivity(...)

            /**
                \brief Can be used to form the matrix of Lagrange multipliers (for debugging purposes).

                \note In the main ObjectiveSensitivity(...) function the real residual might be used. Here
                only the residual based on the factorization is used.
            */
            void ObjectiveSensitivity(Index ObjIndex)
            {
                Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
                Index &ColDim = FirstColIndex; // ColDim will be used to indicate column dimension (I use it for clarity)

                // Even though this computation can be improved, it is very cheap and it is convenient to perform here
                Index nLambda = 0; // Number of constraints that may influence objective ObjIndex (excluding "variable fixing" constraints)
                Index nRank   = 0; // Total rank of objectives before objective ObjIndex
                for (Index k=0; k<ObjIndex; k++)
                {
                    nLambda += obj_info[k].dim;
                    nRank   += obj_info[k].rank;
                }
                nLambda += obj_info[ObjIndex].dim;

                // ---------------------------------------------------------------------------------
                dWorkspace.head(nVarFixed + nLambda + nRank + nVarFixed).setZero();
                dVectorBlockType LambdaFixed(dWorkspace,                   0, nVarFixed);
                dVectorBlockType      Lambda(dWorkspace,           nVarFixed, nLambda);
                dVectorBlockType         rhs(dWorkspace, nLambda + nVarFixed, nRank+nVarFixed);
                // ---------------------------------------------------------------------------------

                FirstRowIndex = obj_info[ObjIndex].first_row_index;
                FirstColIndex = obj_info[ObjIndex].first_col_index;
                ObjDim        = obj_info[ObjIndex].dim;
                ObjRank       = obj_info[ObjIndex].rank;

                bool compute_residual_from_factorization = true;
                if (parameters.regularization_type == 7) // WARNING: TESTING
                {
                    compute_residual_from_factorization = false;
                }

                if (compute_residual_from_factorization)
                {
                    // Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed

                    // copy only what is needed to compute the residual v = A*x-b (i.e., -y_hat)
                    Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = \
                        -LOD.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);

                    // compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
                    Lambda.segment(FirstRowIndex, ObjDim)
                        .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                      FirstColIndex,
                                                                      ObjDim,
                                                                      ObjRank),
                                                            hh_scalars.segment(FirstRowIndex,ObjDim)));
                }
                else
                {
                    initialize_rhs(ObjIndex,rhs);
                    Lambda.segment(FirstRowIndex, ObjDim) = residual_mu.segment(FirstRowIndex,ObjDim);
                }

                if (ObjIndex>0) // the first objective has only Lagrange multipliers equal to the optimal residual
                {
                    // e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}
                    rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                        Lambda.segment(FirstRowIndex, ObjDim);

                    //for (int k=ObjIndex-1; k>=0; k--)
                    for (Index k=ObjIndex; k--; ) //ObjIndex-1, ..., 0 (i.e., for all objectives before ObjIndex)
                    {
                        FirstRowIndex = obj_info[k].first_row_index;
                        FirstColIndex = obj_info[k].first_col_index;
                        ObjDim        = obj_info[k].dim;
                        ObjRank       = obj_info[k].rank;

                        // Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed

                        Lambda.segment(FirstRowIndex, ObjRank) = rhs.segment(FirstColIndex, ObjRank);

                        // apply Q_k' on the left
                        Lambda.segment(FirstRowIndex, ObjDim)
                            .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                          FirstColIndex,
                                                                          ObjDim,
                                                                          ObjRank),
                                                                hh_scalars.segment(FirstRowIndex,ObjDim)));

                        rhs.head(ColDim).noalias() -= LOD.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                            Lambda.segment(FirstRowIndex, ObjDim);

                    } // END for (Index k=ObjIndex-1; k>=0; k--)
                }

                if (nVarFixed>0) // Handle fixed variables (if any)
                {
                    LambdaFixed = -LOD.block(0, 0, nLambda, nVarFixed).transpose() * Lambda;
                }

                rhs.setZero(); // for convenience (easier to analyze the Lagrange multipliers by hand)

            } // END ObjectiveSensitivity(...)


            /**
               \brief Output all ctr whose Lambdas have wrong sign
            */
            void findDescentDirection(Index ObjIndex,
                                      int FirstRowIndex,
                                      Index ObjDim,
                                      const dVectorType& lambda,
                                      RealScalar tol_wrong_sign_lambda,
                                      RealScalar tol_correct_sign_lambda,
                                      std::vector<ConstraintInfo> &ctr_wrong_sign)
            {
                RealScalar aLambda;
                Index ind;
                ConstraintActivationType *aCtrType; // aCtrType cannot be CTR_INACTIVE

                for (Index k=0; k<ObjDim; k++)
                {
                    if (FirstRowIndex < 0) // handle fixed variables
                    {
                        ind      = k;
                        aCtrType = &fixed_var_type[ind];
                    }
                    else                   // handle general constraints
                    {
                        ind      = FirstRowIndex+k;
                        aCtrType = &ctr_type[ind];
                    }

                    if (*aCtrType != CTR_ACTIVE_EQ && *aCtrType != CORRECT_SIGN_OF_LAMBDA)
                    {
                        aLambda = lambda.coeff(ind);

                        if (*aCtrType == CTR_ACTIVE_LB)
                        {
                            aLambda = -aLambda;
                        }

                        if (aLambda > tol_correct_sign_lambda)
                        {
                            *aCtrType = CORRECT_SIGN_OF_LAMBDA;
                        }
                        else if (aLambda < -tol_wrong_sign_lambda)
                        {
                            ctr_wrong_sign.push_back(ConstraintInfo((int)ObjIndex,(int)k));
                        }
                    }
                }
            }


            /**
               \brief Given a vector of Lagrange multipliers, determine the largest (in absolute value)
               multiplier with a wrong sign.

               \param[in]     FirstRowIndex  Index of first element of the current objective (if FirstRowIndex < 0 handle fixed variables)
               \param[in]     ObjDim         Number of constraints.
               \param[in,out] maxAbsValue    Largest (in absolute value) multiplier with wrong sign
               \param[in,out] CtrIndex       Index of largest (in absolute value) multiplier with wrong sign.
               \param[in]     lambda         Vector of lagrange multipliers.
               \param[in]     tol_wrong_sign_lambda   Absolute value of Lagrange multiplier to be considered with "wrong" sign.
               \param[in]     tol_correct_sign_lambda Absolute value of Lagrange multiplier to be considered with "correct" sign.

               \note Lagrange multipliers in the interval (-tol_wrong_sign_lambda --- 0 --- tol_correct_sign_lambda) are considered equal to zero.

               \note Using tol_correct_sign_lambda = 0 is not a good idea in general. This might lead to
               situations where in order to decrease the norm of the residual of a higher-level task
               with e.g., 1e-12, the solver might worsen the nor of the residual of lower-level taks
               with a lot more (e.g., 10).

               \return true if there are multipliers with a wrong sign whose absolute value is larger
               than the largest multiplier with a wrong sign from previous groups of Lagrange
               multipliers.
            */
            bool findDescentDirection(int FirstRowIndex,
                                      Index ObjDim,
                                      RealScalar &maxAbsValue,   // modified
                                      Index &CtrIndex,           // modified
                                      const dVectorType& lambda,
                                      RealScalar tol_wrong_sign_lambda,
                                      RealScalar tol_correct_sign_lambda)
            {
                bool FoundBetterDescentDirection = false;
                RealScalar aLambda;
                Index ind;
                ConstraintActivationType *aCtrType; // aCtrType cannot be CTR_INACTIVE

                for (Index k=0; k<ObjDim; k++)
                {
                    if (FirstRowIndex < 0) // handle fixed variables
                    {
                        ind      = k;
                        aCtrType = &fixed_var_type[ind];
                    }
                    else                   // handle general constraints
                    {
                        ind      = FirstRowIndex+k;
                        aCtrType = &ctr_type[ind];
                    }

                    if (*aCtrType != CTR_ACTIVE_EQ && *aCtrType != CORRECT_SIGN_OF_LAMBDA)
                    {
                        aLambda = lambda.coeff(ind);

                        if (*aCtrType == CTR_ACTIVE_LB)
                        {
                            aLambda = -aLambda;
                        }

                        if (aLambda > tol_correct_sign_lambda)
                        {
                            *aCtrType = CORRECT_SIGN_OF_LAMBDA;
                        }
                        else if (aLambda < -tol_wrong_sign_lambda)
                        {
                            if (aLambda < maxAbsValue) // heuristics: find the multiplier with largest absolute value
                            {
                                FoundBetterDescentDirection = true;
                                maxAbsValue = aLambda;
                                CtrIndex    = k;
                            }
                        }
                    }
                }

                return FoundBetterDescentDirection;
            }

            /**
               \brief Back-solve accounting for the zero blocks due to singular constraints

               \note Here I use x as a temporary variable so that at each level of the recursion I can
               perform only one matrix-vector product. The overhead of copying the RHS to x is smaller
               compared to the gain from performing only one matrix-vector product.

               \verbatim
               | A  A2  A3 | |x1|   |b1|
               | 0  B   B3 |*|x2| = |b2|
               | 0  0   C  | |x3|   |b3|

               x3 = C\b3

               x2 = b2 - B3*x3
               x2 = B\x2;

               x1 = b1 - [A2,A3]*[x2;x3]
               x1 = A\x1;
               \endverbatim

               \note The current implementation does not handle successive objectives whose constraints
               are not singular in an optimal way (in which case directly Eigen's triangular solver
               could be used). As a result the code is much simpler, and just a bit slower (~1% for the
               problems I have tested) from the built-in solvers.
            */
            void solve()
            {
                Index ObjRank, AccumulatedRanks = 0;
                //for(Index k=nObj-1; k>=0; k--)
                for (Index k=nObj; k--; ) //nObj-1, ..., 0.
                {
                    ObjRank = obj_info[k].rank;
                    if (ObjRank > 0)
                    {
                        dVectorBlockType x_k(x, obj_info[k].first_col_index, ObjRank);

                        x_k = LOD.col(nVar).segment(obj_info[k].first_row_index, ObjRank);

                        if (AccumulatedRanks > 0) // Do not enter here the first time ObjRank != 0
                        {
                            x_k.noalias() -= LOD.block(obj_info[k].first_row_index,
                                                       obj_info[k+1].first_col_index,
                                                       ObjRank,
                                                       AccumulatedRanks) * x.segment(obj_info[k+1].first_col_index, AccumulatedRanks);
                        }

                        LOD.block(obj_info[k].first_row_index, obj_info[k].first_col_index, ObjRank, ObjRank)
                            .triangularView<Eigen::Upper>().solveInPlace(x_k);

                        AccumulatedRanks += ObjRank;
                    }
                }

                // Apply permutation
                x = P*x;
            }

            /**
               \brief Compute the least-norm solution.

               \note using Givens rotations
            */
            void solveLeastNorm_1()
            {
                Index FirstRowIndex, FirstColIndex, ObjRank, counter = 0;
                Index nVarRank = 0; // number of variables determined by rank([R,T])

                // -------------------------------------------------------------------------
                // determine dimensions
                // -------------------------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nVarRank += obj_info[ObjIndex].rank;
                }

                Index nVarFree = nVar - (nVarRank + nVarFixed);

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType  R(array, 0,        0, nVarRank, nVarRank);
                dBlockType  T(array, 0, nVarRank, nVarRank, nVarFree);
                dBlockType RT(array, 0,        0, nVarRank, nVarRank+nVarFree); // for convenience

                dBlockType2Vector rhs(array, 0, nVar, nVarRank + nVarFree, 1);
                rhs.tail(nVarFree).setZero(); // important

                // -------------------------------------------------------------------------
                // copy stuff
                // -------------------------------------------------------------------------
                Index col_dim = nVarRank + nVarFree;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    RT.block(counter, counter, ObjRank, col_dim)
                        .triangularView<Eigen::Upper>() = LOD.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim);

                    rhs.segment(counter,ObjRank) = LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                    counter += ObjRank;
                    col_dim -= ObjRank;
                }

                // -------------------------------------------------------------------------
                // zero-out the redundant part (by applying Givens rotations on the right)
                // -------------------------------------------------------------------------
                GivensRotationSequence gs(nVarFree*nVarRank);
                for (Index i=0; i<nVarFree; i++)
                {
                    //for (int j=nVarRank-1; j>=0; j--)
                    for (Index j=nVarRank; j--; ) //nVarRank-1, ..., 0.
                    {
                        GivensRotation GR(RT.coeffRef(j,j),RT.coeffRef(j,nVarRank+i),j,nVarRank+i);
                        RT.topRows(j+1).applyOnTheRight(GR.i,GR.j,GR.G);

                        gs.push(GR);
                    }
                }

                // -------------------------------------------------------------------------
                // backward substitution
                // -------------------------------------------------------------------------
                R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(rhs.head(nVarRank));

                // -------------------------------------------------------------------------
                // apply sequence of Givens rotations on the RHS vector
                // -------------------------------------------------------------------------
                //for (Index i=gs.size()-1; i>=0; i--)
                for (Index i=gs.size(); i--; ) //gs.size()-1, ..., 0.
                {
                    rhs.applyOnTheLeft(gs.get_i(i), gs.get_j(i), gs.get(i));
                }

                // -------------------------------------------------------------------------
                // Apply permutation
                // -------------------------------------------------------------------------
                x.tail(nVarRank + nVarFree) = rhs; // x.head(nVarFixed) contain the fixed variables
                x = P*x;
            }

            /**
               \brief Compute the least-norm solution.

               \note using the normal equations
            */
            void solveLeastNorm_2()
            {
                Index FirstRowIndex, FirstColIndex, ObjRank, counter = 0;
                Index nVarRank = 0; // number of variables determined by rank([R,T])

                // -------------------------------------------------------------------------
                // determine dimensions
                // -------------------------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nVarRank += obj_info[ObjIndex].rank;
                }

                Index nVarFree = nVar - (nVarRank + nVarFixed);

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType  R(array, 0,        0, nVarRank, nVarRank);
                dBlockType  T(array, 0, nVarRank, nVarRank, nVarFree+1);
                dBlockType RT(array, 0,        0, nVarRank, nVarRank+nVarFree+1); // for convenience

                dBlockType  D(array, nVarRank, 0, nVarFree, nVarFree);
                dBlockType2Vector d(array, 0, nVar, nVarFree, 1);

                // -------------------------------------------------------------------------
                // copy stuff
                // -------------------------------------------------------------------------
                Index col_dim = nVarRank + nVarFree;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    RT.block(counter, counter, ObjRank, col_dim+1)
                        .triangularView<Eigen::Upper>() = LOD.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim+1);

                    counter += ObjRank;
                    col_dim -= ObjRank;
                }

                // no need to negate
                R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(T);

                D.triangularView<Eigen::Lower>() = T.leftCols(nVarFree).transpose()*T.leftCols(nVarFree);
                for (Index i=0; i<nVarFree; i++)
                {
                    D.coeffRef(i,i) += 1.0;
                }

                d = T.leftCols(nVarFree).transpose() * T.col(nVarFree);

                Eigen::LLT<dMatrixType> chol(D);
                x.tail(nVarFree) = chol.solve(d);

                counter = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    x.segment(nVarFixed+counter,ObjRank) = LOD.col(nVar).segment(FirstRowIndex,ObjRank) -
                        LOD.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);

                    counter += ObjRank;
                }
                R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(x.segment(nVarFixed,nVarRank));

                // -------------------------------------------------------------------------
                // Apply permutation
                // -------------------------------------------------------------------------
                x = P*x;
            }

            /**
               \brief Compute the least-norm solution.

               \note using the normal equations (and reusing the basis constructed for Tikhonov
               regularization). Hence in order to use this, one has to set regularization_type =
               REGULARIZATION_TIKHONOV and obj_info[:].regularization_factor = 0;
            */
            void solveLeastNorm_3()
            {
                Index FirstRowIndex, ObjRank, counter = 0;
                Index nVarRank = 0; // number of variables determined by rank([R,T])

                // -------------------------------------------------------------------------
                // determine dimensions
                // -------------------------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nVarRank += obj_info[ObjIndex].rank;
                }

                Index nVarFree = nVar - (nVarRank + nVarFixed);

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType iR(null_space, 0,          nVarFixed, nVarRank,   nVarRank); // inv(R)
                dBlockType  T(null_space, 0, nVarFixed+nVarRank, nVarRank, nVarFree+1); // inv(R)*T
                dBlockType  D(     array, 0,                  0, nVarFree,   nVarFree);

                dBlockType2Vector d(array, 0, nVar, nVarFree, 1);

                // -------------------------------------------------------------------------

                D.triangularView<Eigen::Lower>() = T.leftCols(nVarFree).transpose()*T.leftCols(nVarFree);
                for (Index i=0; i<nVarFree; i++)
                {
                    D.coeffRef(i,i) += 1.0;
                }

                d = T.leftCols(nVarFree).transpose() * T.col(nVarFree);

                Eigen::LLT<dMatrixType> chol(D);
                x.tail(nVarFree) = chol.solve(d);

                counter = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    x.segment(nVarFixed+counter,ObjRank) = LOD.col(nVar).segment(FirstRowIndex,ObjRank) -
                        LOD.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);

                    counter += ObjRank;
                }
                x.segment(nVarFixed,nVarRank) = iR.triangularView<Eigen::Upper>() * x.segment(nVarFixed,nVarRank);

                // -------------------------------------------------------------------------
                // Apply permutation
                // -------------------------------------------------------------------------
                x = P*x;
            }


            /**
               \brief Compute the solution that minimizes ||M.leftCols(nVar)*x - M.col(nVar)||^2.

               \note using the normal equations

               \note M is copied (because I modify it below)
            */
            void solveGeneralNorm(dMatrixType M)
            {
                Index FirstRowIndex, FirstColIndex, ObjRank, counter = 0;
                Index nVarRank = 0; // number of variables determined by rank([R,T])

                M = M*P; // permute columns

                // -------------------------------------------------------------------------
                // determine dimensions
                // -------------------------------------------------------------------------
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nVarRank += obj_info[ObjIndex].rank;
                }

                Index nVarFree = nVar - (nVarRank + nVarFixed);

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType  R(array, 0,        0, nVarRank, nVarRank);
                dBlockType  T(array, 0, nVarRank, nVarRank, nVarFree+1);
                dBlockType RT(array, 0,        0, nVarRank, nVarRank+nVarFree+1); // for convenience

                dBlockType  D(array, nVarRank, 0, nVarFree, nVarFree);
                dBlockType2Vector d(array, 0, nVar, nVarFree, 1);

                dBlockType LB(M, 0,        0, nVar, nVarRank);
                dBlockType TB(M, 0, nVarRank, nVar, nVarFree+1);

                // -------------------------------------------------------------------------
                // copy stuff
                // -------------------------------------------------------------------------
                Index col_dim = nVarRank + nVarFree;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    RT.block(counter, counter, ObjRank, col_dim+1)
                        .triangularView<Eigen::Upper>() = LOD.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim+1);

                    counter += ObjRank;
                    col_dim -= ObjRank;
                }

                // -------------------------------------------------------------------------

                R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(LB);
                TB.noalias() -= LB * T;

                D.triangularView<Eigen::Lower>() = TB.leftCols(nVarFree).transpose()*TB.leftCols(nVarFree);
                d = TB.leftCols(nVarFree).transpose() * TB.col(nVarFree);

                Eigen::LLT<dMatrixType> chol(D);
                x.tail(nVarFree) = chol.solve(d);

                counter = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjRank       = obj_info[ObjIndex].rank;

                    x.segment(nVarFixed+counter,ObjRank) = LOD.col(nVar).segment(FirstRowIndex,ObjRank) -
                        LOD.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);

                    counter += ObjRank;
                }
                R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(x.segment(nVarFixed,nVarRank));

                // -------------------------------------------------------------------------
                // Apply permutation
                // -------------------------------------------------------------------------
                x = P*x;
            }

            // =================================================================================================
            // set & get
            // =================================================================================================

            /**
               \brief Declare a variable as fixed

               \param[in] VarIndex Index of variable: x(VarIndex) will be fixed
               \param[in] VarValue Fix x(VarIndex) = VarValue
               \param[in] type     CTR_ACTIVE_LB if x(VarIndex) is fixed at its lower bound, CTR_ACTIVE_UB if
               x(VarIndex) is fixed at its upper bound in LexLSI

               \note When we don't care about the type (i.e., we will not use the associated Lagrange
               multipliers), CTR_ACTIVE_UB is (arbitrary) chosen as a default value for type.

            */
            void fixVariable(Index VarIndex, RealScalar VarValue, ConstraintActivationType type = CTR_ACTIVE_UB)
            {
                fixed_var_index(nVarFixedInit) = VarIndex;
                x(nVarFixedInit)               = VarValue;
                fixed_var_type[nVarFixedInit]  = type;

                nVarFixedInit++;
            }

            /**
                \brief Declare some of the variables as fixed

                \param[in] nVarFixed_ Number of fixed variables
                \param[in] VarIndex   Indexes of variables to fix
                \param[in] VarValue   Values of the fixed variables
                \param[in] type       can be CTR_ACTIVE_LB or CTR_ACTIVE_UB
            */
            void fixVariables(Index nVarFixed_, Index *VarIndex, RealScalar *VarValue, ConstraintActivationType *type)
            {
                if (nVarFixed_ > nVar)
                {
                    throw Exception("Cannot fix more than nVar variables");
                }
                else
                {
                    nVarFixed = nVarFixed_;
                }

                fixed_var_index.resize(nVarFixed);
                fixed_var_type.resize(nVarFixed);

                // copy data
                fixed_var_index     = Eigen::Map<iVectorType>(VarIndex, nVarFixed);
                x.head(nVarFixed) = Eigen::Map<dVectorType>(VarValue, nVarFixed); // x will be permuted later
                fixed_var_type.assign(type, type+nVarFixed);
            }

            /**
                \brief Set dimension of objectives

                \param[in] ObjDim_ Number of constraints involved in each objective
            */
            void setObjDim(Index *ObjDim_)
            {
                nCtr = 0;
                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    nCtr += ObjDim_[ObjIndex];

                    obj_info[ObjIndex].dim = ObjDim_[ObjIndex];
                    if (ObjIndex > 0)
                    {
                        obj_info[ObjIndex].first_row_index = obj_info[ObjIndex-1].first_row_index + obj_info[ObjIndex-1].dim;
                    }
                }

                initialize();
            }


            /**
                \brief Set number of fixed variables

                \param[in] nVarFixed_ Number of fixed variables
            */
            void setFixedVariablesCount(Index nVarFixed_)
            {
                if (nVarFixed_ > nVar)
                {
                    throw Exception("Cannot fix more than nVar variables");
                }
                else
                {
                    nVarFixed = nVarFixed_;
                }

                fixed_var_index.resize(nVarFixed);
                fixed_var_type.resize(nVarFixed);
            }


            /**
               \brief Sets parameters
            */
            void setParameters(const ParametersLexLSE &parameters_)
            {
                parameters = parameters_;
            }

            /**
                \brief Set (a non-negative) regularization factor for objective ObjIndex

                \note Note that fixed variables are not counted as an objective
            */
            void setRegularizationFactor(Index ObjIndex, RealScalar factor)
            {
                /// @todo: check whether ObjIndex and factor make sense.

                obj_info[ObjIndex].regularization_factor = factor;
            }

            /**
                \brief Get number of fixed variables
            */
            Index getFixedVariablesCount()
            {
                return nVarFixed;
            }

            /**
                \brief Get indexes of fixed variables
            */
            iVectorType& getFixedVarIndex()
            {
                return fixed_var_index;
            }

            /**
               \brief Return #TotalRank (should be called after factorizing)
            */
            Index getTotalRank()
            {
                return TotalRank;
            }

            /**
                \brief set a random problem [A,RHS]
            */
            void setProblem(const dMatrixType& data)
            {
                LOD = data;
            }

            /**
                \brief Set data of objective ObjIndex [A,RHS] - the data is copied

                \param[in] ObjIndex Index of objective
                \param[in] data     data (including LHS & RHS)
            */
            void setData(Index ObjIndex, const dMatrixType& data)
            {
                if (ObjIndex >= nObj)
                {
                    throw Exception("ObjIndex >= nObj");
                }

                LOD.block(obj_info[ObjIndex].first_row_index,0,obj_info[ObjIndex].dim,nVar+1) = data;
            }

            /**
                \brief Set one constraint

                \param[in] CtrIndex Index of row in LOD (regardless of objective)
                \param[in] row      Constraint vector (LHS)
                \param[in] rhs      RHS vector
            */
            void setCtr(Index CtrIndex, const dRowVectorType& row, RealScalar rhs)
            {
                LOD.row(CtrIndex).head(nVar) = row;
                LOD.coeffRef(CtrIndex,nVar)  = rhs;
            }

            /**
                \brief Set the type of a constraint with index CtrIndex in LexLSE objective ObjIndex
            */
            void setCtrType(Index ObjIndex, Index CtrIndex, ConstraintActivationType type)
            {
                Index FirstRowIndex = obj_info[ObjIndex].first_row_index;
                ctr_type[FirstRowIndex + CtrIndex] = type;
            }

            /**
               \brief Form the residuals (A*x-RHS) through the LOD. The residual of the
               fixed variables is always zero (and is not included).

               \note This function could compute an incorect residual depending on parameters.tol_linear_dependence.
            */
            dVectorType& get_v()
            {
                Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
                dVectorBlockType v(dWorkspace,0,nCtr);

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    FirstRowIndex = obj_info[ObjIndex].first_row_index;
                    FirstColIndex = obj_info[ObjIndex].first_col_index;
                    ObjDim        = obj_info[ObjIndex].dim;
                    ObjRank       = obj_info[ObjIndex].rank;

                    v.segment(FirstRowIndex, ObjRank).setZero(); // Zero-out first ObjRank elements
                    v.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = -LOD.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);

                    v.segment(FirstRowIndex, ObjDim)
                        .applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                      FirstColIndex,
                                                                      ObjDim,
                                                                      ObjRank),
                                                            hh_scalars.segment(FirstRowIndex,ObjDim)));
                }

                return dWorkspace; // return directly a reference to the working array (use dWorkspace.head(nCtr) outside)
            }

            /**
                \brief Return the (primal) solution vector
            */
            dVectorType& get_x()
            {
                return x;
            }

            /**
                \brief Return the number of equations in objective ObjIndex
            */
            Index getDim(Index ObjIndex) const
            {
                return obj_info[ObjIndex].dim;
            }

            /**
                \brief Return the rank of the equations in objective ObjIndex
            */
            Index getRank(Index ObjIndex) const
            {
                return obj_info[ObjIndex].rank;
            }

            Index get_nObj()
            {
                return nObj;
            }

            Index get_nVar()
            {
                return nVar;
            }

            /**
                \brief Return dWorkspace
            */
            dVectorType& getWorkspace()
            {
                return dWorkspace;
            }

            dMatrixType get_lexqr()
            {
                return LOD;
            }

            dMatrixType get_data()
            {
                return PROBLEM_DATA;
            }

            dMatrixType get_X_mu()
            {
                return X_mu;
            }

            dMatrixType get_X_mu_rhs()
            {
                return X_mu_rhs;
            }

            dVectorType get_residual_mu()
            {
                return residual_mu;
            }

            /**
                \brief Reset the factorization (initialize A separately)
            */
            void reset()
            {
                initialize();
                x.head(nVarFixed).setZero();
            }

            // =================================================================================================
            // Utilities
            // =================================================================================================

        private:

            /**
                \brief Initialize some of the fields

                \note Some of the fields have already been initialized in the constructors but it is
                necessary to initialize them again here. This is necessary when an instance of the class
                LexLSE is used to solve multiplie problems (as it is done in LexLSI).
            */
            void initialize()
            {
                nVarFixedInit = 0;
                TotalRank     = 0;

                for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                {
                    obj_info[ObjIndex].rank = 0;
                    obj_info[ObjIndex].first_col_index = 0;
                }

                hh_scalars.setZero();

                null_space.setZero();
                array.setZero();
                X_mu.setZero();
                X_mu_rhs.setZero();
                residual_mu.setZero();

                P.setIdentity();
                x.tail(nVar-nVarFixed).setZero(); // x.head(nVarFixed) has already been initialized in fixVariable(...)
            }

            /**
                \brief Tikhonov regularization (using the normal equations inv(A'*A+I)*A'*b)

                \note fast when column-dimension is small
            */
            void regularize_tikhonov_1(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(LOD,
                              FirstRowIndex, FirstColIndex,
                              ObjRank, ObjRank);

                dBlockType Tk(LOD,
                              FirstRowIndex, FirstColIndex+ObjRank,
                              ObjRank, RemainingColumns);

                // Matrix S_{k-1} (maybe change the name "up")
                dBlockType up(null_space,
                              0, FirstColIndex,
                              FirstColIndex-nVarFixed, RemainingColumns+ObjRank);

                dBlockType  D(array,
                              0, 0,
                              RemainingColumns+ObjRank, RemainingColumns+ObjRank);

                dBlockType2Vector d(array,
                                    0, nVar,
                                    RemainingColumns+ObjRank, 1);

                // -------------------------------------------------------------------------

                // Rk'*Rk
                D.block(0,0,ObjRank,ObjRank)
                    .triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();

                // Tk'*Tk
                D.block(ObjRank,ObjRank,RemainingColumns,RemainingColumns)
                    .triangularView<Eigen::Lower>() = Tk.transpose()*Tk;

                // Tk'*Rk
                D.block(ObjRank,0,RemainingColumns,ObjRank).noalias() = Tk.transpose()*Rk.triangularView<Eigen::Upper>();

                // [Rk, Tk]'*[Rk, Tk] + S_{k-1}'*S_{k-1}
                D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;

                // ... + mu*I
                for (Index i=0; i<RemainingColumns+ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                // ==============================================================================================
                // Rk'*xk_star
                d.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose() * \
                    LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                // Tk'*xk_star
                d.tail(RemainingColumns).noalias() = Tk.transpose() *\
                    LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                // S_{k-1}'*s_{k-1}
                d.noalias() += mu * up.transpose() * null_space.col(nVar).head(FirstColIndex-nVarFixed);

                // ==============================================================================================

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);

                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>() * d.head(ObjRank);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk * d.tail(RemainingColumns);
            }


            /**
                \brief like regularize_tikhonov_1 but for testing stuff

                \todo IMPORTANT: if we don't enter here at iteration i (because the i-th objective
                doesn't have a single pivot), we have to set X_mu.col(i) = X_mu.col(i-1).

                \todo If I endup passig ObjIndex as input argument, there is no need to pass
                FirstRowIndex, FirstColIndex and ObjRank
            */
            void regularize_tikhonov_1_test(Index FirstRowIndex,
                                            Index FirstColIndex,
                                            Index ObjRank,
                                            Index RemainingColumns,
                                            Index ObjIndex)
            {
                Index ObjDim = obj_info[ObjIndex].dim;

                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(LOD,
                              FirstRowIndex, FirstColIndex,
                              ObjRank, ObjRank);

                dBlockType Tk(LOD,
                              FirstRowIndex, FirstColIndex+ObjRank,
                              ObjRank, RemainingColumns);

                // Matrix S_{k-1} (maybe change the name "up")
                dBlockType up(null_space,
                              0, FirstColIndex,
                              FirstColIndex-nVarFixed, RemainingColumns+ObjRank);

                dBlockType  D(array,
                              0, 0,
                              RemainingColumns+ObjRank, RemainingColumns+ObjRank);

                dBlockType2Vector d(array,
                                    0, nVar,
                                    RemainingColumns+ObjRank, 1);

                dBlockType2Vector rhs(LOD, FirstRowIndex, nVar, ObjDim, 1);

                // This block overlaps with the portion of dWorkspacededicated to hhWorkspace, but
                // we can use it here
                dVectorBlockType residual_workspace(dWorkspace, nVar, ObjDim);

                // -------------------------------------------------------------------------

                // Rk'*Rk
                D.block(0,0,ObjRank,ObjRank)
                    .triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();

                // Tk'*Tk
                D.block(ObjRank,ObjRank,RemainingColumns,RemainingColumns)
                    .triangularView<Eigen::Lower>() = Tk.transpose()*Tk;

                // Tk'*Rk
                D.block(ObjRank,0,RemainingColumns,ObjRank).noalias() = Tk.transpose()*Rk.triangularView<Eigen::Upper>();

                // [Rk, Tk]'*[Rk, Tk] + S_{k-1}'*S_{k-1}
                D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;

                // ... + mu*I
                for (Index i=0; i<RemainingColumns+ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                // ==============================================================================================
                // Rk'*xk_star
                d.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose() * rhs.head(ObjRank);

                // Tk'*xk_star
                d.tail(RemainingColumns).noalias() = Tk.transpose() * rhs.head(ObjRank);

                // S_{k-1}'*s_{k-1}
                d.noalias() += mu * up.transpose() * null_space.col(nVar).head(FirstColIndex-nVarFixed);

                // ==============================================================================================

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);
                rhs.head(ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>() * d.head(ObjRank);
                rhs.head(ObjRank).noalias() += Tk * d.tail(RemainingColumns);

                // ==============================================================================================
                // compute residual Q1*[R T]*x - b
                // ==============================================================================================
                residual_workspace.head(ObjRank) = rhs.head(ObjRank);
                residual_workspace.tail(ObjDim-ObjRank).setZero();    // set y_hat = 0
                residual_workspace.applyOnTheLeft(householderSequence(LOD.block(FirstRowIndex,
                                                                                FirstColIndex,
                                                                                ObjDim,
                                                                                ObjRank),
                                                                      hh_scalars.segment(FirstRowIndex,ObjDim)));

                residual_mu.segment(FirstRowIndex,ObjDim) = residual_workspace - residual_mu.segment(FirstRowIndex,ObjDim);
                // ==============================================================================================

                X_mu.col(ObjIndex).tail(RemainingColumns+ObjRank) = d;
                get_intermediate_x(ObjIndex, RemainingColumns+ObjRank);

                // ==============================================================================================
                // permute X_mu.col(ObjIndex)
                // ==============================================================================================
                Index AccumulatedRanks = 0;
                for (Index k=0; k<=ObjIndex; k++)
                {
                    AccumulatedRanks += obj_info[k].rank;
                }

                // equivalent to the code with P_tmp below
                for (Index k=AccumulatedRanks; k--; )
                {
                    Index j = column_permutations.coeff(k);
                    std::swap(X_mu.coeffRef(k,ObjIndex), X_mu.coeffRef(j,ObjIndex));
                }

                /*
                Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P_tmp;
                P_tmp.setIdentity(nVar);
                for (Index k=0; k<AccumulatedRanks; k++)
                {
                    P_tmp.applyTranspositionOnTheRight(k, column_permutations.coeff(k));
                }
                X_mu.col(ObjIndex) = P_tmp * X_mu.col(ObjIndex);
                */
                // ==============================================================================================
            }


            /*
              \brief Initialize right-hand-size in case of regularization

              \verbatim
              |R11'    0    0  | |x1|   |b1|
              |R12'  R22'   0  |*|x2| = |b2|
              |R13'  R23'  R33'| |x3|   |b3|

              Note that |R12 R13| is continuous in memory but |R13' R23'| is not

              ----------------------------
              The basic flow is:
              ----------------------------
              k = 1
              x1 = inv(R11')*b1

              k = 2
              |b2|   |b2|   |R12'|
              |  | = |  | - |    |*x1
              |b3|   |b3|   |R13'|
              x2 = inv(R22')*b2

              k3 = 3
              b3 = b3 - R23'*x2
              x3 = inv(R33')*b3
              ----------------------------

              Depending on ObjIndex, we might not need to update both (b2,b3) for k=2. This is
              related to the fact that we are solving a very particular system (to explain in the
              note).
              \endverbatim
            */
            template <typename Derived>
            void initialize_rhs(Index ObjIndex, Eigen::MatrixBase<Derived> & rhs)
            {
                X_mu_rhs.col(ObjIndex) = X_mu.col(ObjIndex); // copy (I want to preserve X_mu)

                dBlockType2Vector X_mu_k(X_mu_rhs, 0, ObjIndex, nVar, 1);

                aRegularizationFactor = obj_info[ObjIndex].regularization_factor;

                // permute and scale elements (actually, later I don't use the whole vector X_mu_k
                // when computing the Lagrange multipliers, but all this is very cheap ...)
                X_mu_k = P.transpose()*X_mu_k;
                X_mu_k *= -aRegularizationFactor*aRegularizationFactor;

                // last index of interes for objective ObjIndex
                Index last_col_index = obj_info[ObjIndex].first_col_index + obj_info[ObjIndex].rank - 1;
                Index remain_col;

                for (Index k=0; k<=ObjIndex; k++)
                {
                    if (k>0)
                    {
                        // number of remaining columns
                        remain_col = last_col_index - obj_info[k].first_col_index + 1;

                        dBlockType Rkj(LOD,
                                       obj_info[k-1].first_row_index, obj_info[k].first_col_index,
                                       obj_info[k-1].rank, remain_col);

                        X_mu_k.segment(obj_info[k].first_col_index,remain_col).noalias() -= \
                            Rkj.transpose() * X_mu_k.segment(obj_info[k-1].first_col_index,
                                                             obj_info[k-1].rank);
                    }

                    LOD.block(obj_info[k].first_row_index, obj_info[k].first_col_index,
                              obj_info[k].rank, obj_info[k].rank).transpose()
                        .triangularView<Eigen::Lower>()
                        .solveInPlace(X_mu_k.segment(obj_info[k].first_col_index,obj_info[k].rank));
                }
                rhs = X_mu_k.head(rhs.size());
            }
            /*
            template <typename Derived>
            void initialize_rhs(Index ObjIndex, Eigen::MatrixBase<Derived> & rhs)
            {
                dVectorBlockType X_mu_k(dWorkspace, nVar, nVar);
                X_mu_k = X_mu.col(ObjIndex); // copy (I want to preserve X_mu)

                aRegularizationFactor = obj_info[ObjIndex].regularization_factor;

                // permute and scale elements (actually, later I don't use the whole vector X_mu_k
                // when computing the Lagrange multipliers, but all this is very cheap ...)
                X_mu_k = P.transpose()*X_mu_k;
                X_mu_k *= -aRegularizationFactor*aRegularizationFactor;

                // last index of interes for objective ObjIndex
                Index last_col_index = obj_info[ObjIndex].first_col_index + obj_info[ObjIndex].rank - 1;
                Index remain_col;

                for (Index k=0; k<=ObjIndex; k++)
                {
                    if (obj_info[k].rank > 0)
                    {
                        if (k>0)
                        {
                            // number of remaining columns
                            remain_col = last_col_index - obj_info[k].first_col_index + 1;

                            dBlockType Rkj(LOD,
                                           obj_info[k-1].first_row_index, obj_info[k].first_col_index,
                                           obj_info[k-1].rank, remain_col);

                            X_mu_k.segment(obj_info[k].first_col_index,remain_col).noalias() -= \
                                Rkj.transpose() * X_mu_k.segment(obj_info[k-1].first_col_index,
                                                                 obj_info[k-1].rank);
                        }

                        LOD.block(obj_info[k].first_row_index, obj_info[k].first_col_index,
                                  obj_info[k].rank, obj_info[k].rank).transpose()
                            .triangularView<Eigen::Lower>()
                            .solveInPlace(X_mu_k.segment(obj_info[k].first_col_index,obj_info[k].rank));
                    }
                }
                rhs = X_mu_k.head(rhs.size());
            }
            */

            /*
              \brief Compute the solution of the ObjIndex-th regularized problem (needed for
              computing the Lagrange multipliers)
            */
            void get_intermediate_x(Index ObjIndex, Index x_tail_size)
            {
                Index ObjRank, FirstColIndex, FirstRowIndex;

                /*
                  \verbatim
                  | A  A2  A3  A4| |x1|   |b1|
                  | 0  B   B3  B4|*|x2| = |b2|
                  | 0  0   C   C4| |x3|   |b3|
                                   |x4|   |b4|

                  At level 3, we compute [x3;x4] and form
                  b1 = b1 - [A3,A4]*[x3;x4];
                  b2 = b2 - [B3,B4]*[x3;x4];
                  \endverbatim
                */
                if (ObjIndex > 0)
                {
                    for (Index i=0; i<ObjIndex; i++)
                    {
                        Index FirstRowIndex = obj_info[i].first_row_index;
                        Index FirstColIndex = obj_info[i].first_col_index;
                        Index ObjRank       = obj_info[i].rank;

                        dBlockType RTi(LOD,
                                       FirstRowIndex, 0,
                                       ObjRank, nVar);

                        X_mu.col(ObjIndex).segment(FirstColIndex, ObjRank) =
                            LOD.col(nVar).segment(FirstRowIndex, ObjRank) - \
                            RTi.rightCols(x_tail_size) * X_mu.col(ObjIndex).tail(x_tail_size);
                    }
                }

                /*
                  Similar to solve() but starting from ObjIndex-1 and using X_mu instead of x
                */
                Index AccumulatedRanks = 0;
                for (Index k=ObjIndex; k--; ) //ObjIndex-1, ..., 0.
                {
                    FirstRowIndex = obj_info[k].first_row_index;
                    FirstColIndex = obj_info[k].first_col_index;
                    ObjRank       = obj_info[k].rank;
                    if (ObjRank > 0)
                    {
                        dBlockType2Vector x_k(X_mu,
                                              FirstColIndex, ObjIndex,
                                              ObjRank, 1);

                        if (AccumulatedRanks > 0) // Do not enter here the first time ObjRank != 0
                        {
                            dBlockType Rkj(LOD,
                                           FirstRowIndex, obj_info[k+1].first_col_index,
                                           ObjRank, AccumulatedRanks);

                            x_k.noalias() -= Rkj * X_mu.col(ObjIndex).segment(obj_info[k+1].first_col_index,
                                                                              AccumulatedRanks);
                        }

                        LOD.block(FirstRowIndex, FirstColIndex,
                                  ObjRank, ObjRank).triangularView<Eigen::Upper>().solveInPlace(x_k);

                        AccumulatedRanks += ObjRank;
                    }
                }
            }

            /**
                \brief Tikhonov regularization (option: A'*inv(A*A'+I)*b)

                \note fast when row-dimension is small
            */
            void regularize_tikhonov_2(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(       LOD, FirstRowIndex,         FirstColIndex,                         ObjRank,                         ObjRank);
                dBlockType Tk(       LOD, FirstRowIndex, FirstColIndex+ObjRank,                         ObjRank,                RemainingColumns);
                dBlockType up(null_space,             0,         FirstColIndex,         FirstColIndex-nVarFixed,      RemainingColumns + ObjRank);
                dBlockType  D(     array,             0,                     0, FirstColIndex-nVarFixed+ObjRank, FirstColIndex-nVarFixed+ObjRank);

                dBlockType2Vector d(array, 0, nVar, FirstColIndex-nVarFixed+ObjRank, 1);
                // -------------------------------------------------------------------------

                // forming D is very expensive
                D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>()  = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
                D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>() += Tk*Tk.transpose();

                D.block(ObjRank,ObjRank,FirstColIndex-nVarFixed,FirstColIndex-nVarFixed).triangularView<Eigen::Lower>() = mu*up*up.transpose();

                D.block(ObjRank, 0, FirstColIndex-nVarFixed, ObjRank).noalias()  =
                    aRegularizationFactor*up.leftCols(ObjRank)*Rk.triangularView<Eigen::Upper>().transpose();

                D.block(ObjRank, 0, FirstColIndex-nVarFixed, ObjRank).noalias() +=
                    aRegularizationFactor*up.rightCols(RemainingColumns)*Tk.transpose();

                for (Index i=0; i<FirstColIndex-nVarFixed+ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                d.head(ObjRank) = LOD.col(nVar).segment(FirstRowIndex,ObjRank);
                d.segment(ObjRank,FirstColIndex-nVarFixed) = aRegularizationFactor*null_space.col(nVar).head(FirstColIndex-nVarFixed);

                // -------------------------------------------------------------------------

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);

                // -------------------------------------------------------------------------

                for (Index i=0; i<FirstColIndex-nVarFixed+ObjRank; i++)
                {
                    D.coeffRef(i,i) -= mu;
                }

                d = D.selfadjointView<Eigen::Lower>() * d;
                LOD.col(nVar).segment(FirstRowIndex,ObjRank) = d.head(ObjRank);
            }

            /**
                \brief Tikhonov regularization (basic variables)
            */
            void regularize_R(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType  Rk(       LOD, FirstRowIndex, FirstColIndex,                 ObjRank, ObjRank);
                dBlockType  up(null_space,             0, FirstColIndex, FirstColIndex-nVarFixed, ObjRank);
                dBlockType   D(     array,             0,             0,                 ObjRank, ObjRank);

                dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

                // -------------------------------------------------------------------------

                D.triangularView<Eigen::Lower>()  = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
                D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;
                for (Index i=0; i<ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                d.noalias()  = up.transpose() * null_space.col(nVar).head(FirstColIndex-nVarFixed);
                d *= mu;
                d.noalias() += Rk.triangularView<Eigen::Upper>().transpose()*LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                // -------------------------------------------------------------------------

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = Rk.triangularView<Eigen::Upper>() * d;
            }

            /**
                \brief Tikhonov regularization (basic variables no Z)
            */
            void regularize_R_NO_Z(Index FirstRowIndex, Index FirstColIndex, Index ObjRank)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk( LOD, FirstRowIndex, FirstColIndex, ObjRank, ObjRank);
                dBlockType  D(array,            0,             0, ObjRank, ObjRank);

                dBlockType2Vector d(array, 0, nVar, ObjRank, 1);
                // -------------------------------------------------------------------------

                D.triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
                for (Index i=0; i<ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                d.noalias() = Rk.triangularView<Eigen::Upper>().transpose()*LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                // -------------------------------------------------------------------------

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = Rk.triangularView<Eigen::Upper>() * d;
            }

            /**
                \brief RT regularization (no Z): [R,T;I]*x = [b;0]
            */
            void regularize_RT_NO_Z(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(  LOD, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
                dBlockType Tk(  LOD, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);
                dBlockType  D(array,             0,                     0, ObjRank,          ObjRank);

                dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

                // -------------------------------------------------------------------------

                D.triangularView<Eigen::Lower>() = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
                D.triangularView<Eigen::Lower>() += Tk*Tk.transpose();

                for (Index i=0; i<ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                d = LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);

                // -------------------------------------------------------------------------

                for (Index i=0; i<ObjRank; i++)
                {
                    D.coeffRef(i,i) -= mu;
                }
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = D.selfadjointView<Eigen::Lower>() * d;
            }

            /**
               \brief For testing purposes
            */
            void regularize_test(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                // just scale the RHS (not a good idea)
                LOD.col(nVar).segment(FirstRowIndex,ObjRank) *= aRegularizationFactor;
            }

            /**
                \brief Tikhonov regularization using CGLS
            */
            void regularize_tikhonov_CG(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(LOD, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
                dBlockType Tk(LOD, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

                // ----------------------------------------------------------------------------------------------
                // generate x0
                // ----------------------------------------------------------------------------------------------
                dVectorType sol_x(ObjRank+RemainingColumns);
                sol_x.setZero();
                // ----------------------------------------------------------------------------------------------

                cg_tikhonov(sol_x, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol_x.tail(RemainingColumns);
            }

            /**
                \brief Tikhonov regularization using CGLS

                \note hot-start from RT_NO_Z
            */
            void regularize_tikhonov_CG_x0(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar mu = aRegularizationFactor*aRegularizationFactor;

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(  LOD, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
                dBlockType Tk(  LOD, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);
                dBlockType  D(array,             0,                     0, ObjRank,          ObjRank);

                dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

                // ----------------------------------------------------------------------------------------------
                // generate x0
                // ----------------------------------------------------------------------------------------------
                D.triangularView<Eigen::Lower>() = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
                D.triangularView<Eigen::Lower>() += Tk*Tk.transpose();

                for (Index i=0; i<ObjRank; i++)
                {
                    D.coeffRef(i,i) += mu;
                }

                d = LOD.col(nVar).segment(FirstRowIndex,ObjRank);

                Eigen::LLT<dMatrixType> chol(D);
                chol.solveInPlace(d);

                dVectorType sol(ObjRank+RemainingColumns);
                sol.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose() * d;
                sol.tail(RemainingColumns) = Tk.transpose() * d;
                // ----------------------------------------------------------------------------------------------

                cg_tikhonov(sol, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol.head(ObjRank);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol.tail(RemainingColumns);
            }

            /**
                \brief RT_NO_Z regularization using CGLS
            */
            void regularize_RT_NO_Z_CG(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(LOD, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
                dBlockType Tk(LOD, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

                // ----------------------------------------------------------------------------------------------
                // generate x0
                // ----------------------------------------------------------------------------------------------
                dVectorType sol(ObjRank+RemainingColumns);
                sol.setZero();
                // ----------------------------------------------------------------------------------------------

                cg_RT(sol, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol.head(ObjRank);
                LOD.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol.tail(RemainingColumns);
            }

            /*
              --------------------------------------------------------------------------
              ObjIndex = k (col_dim: rk + RemainingColumns)
              --------------------------------------------------------------------------
              | Rk Tk | yk | row_dim: rk
              |  Sk   | sk | row_dim: FirstColIndex = ...
              |  Ik   | 0  | row_dim: rk + RemainingColumns
              --------------------------------------------------------------------------
            */
            Index cg_tikhonov(dVectorType &sol_x,
                              Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar alpha, beta, gamma, gamma_previous;

                RealScalar tol = 1e-12; // todo: user input

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------
                dBlockType Rk(       LOD, FirstRowIndex,         FirstColIndex,                 ObjRank,                    ObjRank);
                dBlockType Tk(       LOD, FirstRowIndex, FirstColIndex+ObjRank,                 ObjRank,           RemainingColumns);
                dBlockType Sk(null_space,             0,         FirstColIndex, FirstColIndex-nVarFixed, RemainingColumns + ObjRank);

                dBlockType2Vector r(rqsp_work, 0, 0, ObjRank + FirstColIndex-nVarFixed + ObjRank + RemainingColumns, 1);
                dBlockType2Vector q(rqsp_work, 0, 1, ObjRank + FirstColIndex-nVarFixed + ObjRank + RemainingColumns, 1);
                dBlockType2Vector s(rqsp_work, 0, 2,                                     ObjRank + RemainingColumns, 1);
                dBlockType2Vector p(rqsp_work, 0, 3,                                     ObjRank + RemainingColumns, 1);

                // ------------------------------------------------------------------------------------------------
                /*
                  | yk |   | Rk Tk |
                  r = | sk | - |  Sk   | * x
                  |  0 |   |  Ik   |
                */
                // ------------------------------------------------------------------------------------------------
                r.head(ObjRank).noalias()  = LOD.col(nVar).segment(FirstRowIndex,ObjRank) - Tk*sol_x.tail(RemainingColumns);
                r.head(ObjRank).noalias() -= Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);

                r.segment(ObjRank,FirstColIndex-nVarFixed).noalias() = null_space.col(nVar).head(FirstColIndex-nVarFixed) - Sk * sol_x;
                r.segment(ObjRank,FirstColIndex-nVarFixed)          *= aRegularizationFactor;

                r.tail(ObjRank+RemainingColumns).noalias() = -aRegularizationFactor*sol_x;
                // ------------------------------------------------------------------------------------------------
                /*
                  | Rk'       |   | r1 |
                  s = |     Sk' I | * | r2 |
                  | Tk'       |   | r3 |
                */
                // ------------------------------------------------------------------------------------------------
                s.noalias() = Sk.transpose() * r.segment(ObjRank,FirstColIndex-nVarFixed) + r.tail(ObjRank+RemainingColumns);
                s *= aRegularizationFactor;

                s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
                s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
                // ------------------------------------------------------------------------------------------------
                p = s;
                // ------------------------------------------------------------------------------------------------
                gamma = s.squaredNorm();
                // ------------------------------------------------------------------------------------------------

                Index iter = 0;
                while ( (std::sqrt(gamma) > tol) && (iter < parameters.max_number_of_CG_iterations) )
                {
                    // ------------------------------------------------------------------------------------------------
                    // q = [Rk Tk; Sk ; Ik]*p
                    // ------------------------------------------------------------------------------------------------
                    q.head(ObjRank).noalias()  = Tk*p.tail(RemainingColumns);
                    q.head(ObjRank).noalias() += Rk.triangularView<Eigen::Upper>()*p.head(ObjRank);

                    q.segment(ObjRank,FirstColIndex-nVarFixed).noalias() = Sk * p;
                    q.segment(ObjRank,FirstColIndex-nVarFixed) *= aRegularizationFactor;

                    q.tail(ObjRank+RemainingColumns).noalias() = aRegularizationFactor*p;
                    // ------------------------------------------------------------------------------------------------
                    alpha = gamma/q.squaredNorm();
                    sol_x += alpha*p;
                    r -= alpha*q;
                    // ------------------------------------------------------------------------------------------------
                    // S = [Rk Tk; Sk ; Ik]'*r
                    // ------------------------------------------------------------------------------------------------
                    s.noalias() = Sk.transpose() * r.segment(ObjRank,FirstColIndex-nVarFixed) + r.tail(ObjRank+RemainingColumns);
                    s *= aRegularizationFactor;

                    s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
                    s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
                    // ------------------------------------------------------------------------------------------------
                    gamma_previous = gamma;
                    gamma = s.squaredNorm();
                    beta = gamma/gamma_previous;
                    p = s + beta*p;
                    // ------------------------------------------------------------------------------------------------
                    iter++;
                }

                return iter;
            }

            /*
              --------------------------------------------------------------------------
              ObjIndex = k (col_dim: rk + RemainingColumns)
              --------------------------------------------------------------------------
              | Rk Tk | yk | row_dim: rk
              |  Ik   | 0  | row_dim: rk + RemainingColumns
              --------------------------------------------------------------------------
            */
            Index cg_RT(dVectorType &sol_x,
                        Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                RealScalar alpha, beta, gamma, gamma_previous;

                RealScalar tol = 1e-12; // todo: user input

                // -------------------------------------------------------------------------
                // create blocks
                // -------------------------------------------------------------------------

                dBlockType Rk( LOD, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
                dBlockType Tk( LOD, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

                dBlockType2Vector r(rqsp_work, 0, 0, ObjRank + FirstColIndex + ObjRank + RemainingColumns, 1);
                dBlockType2Vector q(rqsp_work, 0, 1, ObjRank + FirstColIndex + ObjRank + RemainingColumns, 1);
                dBlockType2Vector s(rqsp_work, 0, 2,                           ObjRank + RemainingColumns, 1);
                dBlockType2Vector p(rqsp_work, 0, 3,                           ObjRank + RemainingColumns, 1);

                // ------------------------------------------------------------------------------------------------
                /*
                  | yk |   | Rk Tk |
                  r = |    | - |       | * x
                  |  0 |   |  Ik   |
                */
                // ------------------------------------------------------------------------------------------------
                r.head(ObjRank).noalias()  = LOD.col(nVar).segment(FirstRowIndex,ObjRank) - Tk*sol_x.tail(RemainingColumns);
                r.head(ObjRank).noalias() -= Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);

                r.tail(ObjRank+RemainingColumns).noalias() = -aRegularizationFactor*sol_x;
                // ------------------------------------------------------------------------------------------------
                /*
                  | Rk'   |   | r1 |
                  s = |     I | * | r2 |
                  | Tk'   |
                */
                // ------------------------------------------------------------------------------------------------
                s.noalias() = aRegularizationFactor * r.tail(ObjRank+RemainingColumns);

                s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
                s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
                // ------------------------------------------------------------------------------------------------
                p = s;
                // ------------------------------------------------------------------------------------------------
                gamma = s.squaredNorm();
                // ------------------------------------------------------------------------------------------------

                Index iter = 0;
                while ( (std::sqrt(gamma) > tol) && (iter < parameters.max_number_of_CG_iterations) )
                {
                    // ------------------------------------------------------------------------------------------------
                    // q = [Rk Tk; Ik]*p
                    // ------------------------------------------------------------------------------------------------
                    q.head(ObjRank).noalias()  = Tk*p.tail(RemainingColumns);
                    q.head(ObjRank).noalias() += Rk.triangularView<Eigen::Upper>()*p.head(ObjRank);

                    q.tail(ObjRank+RemainingColumns).noalias() = aRegularizationFactor*p;
                    // ------------------------------------------------------------------------------------------------
                    alpha = gamma/q.squaredNorm();
                    sol_x += alpha*p;
                    r -= alpha*q;
                    // ------------------------------------------------------------------------------------------------
                    // S = [Rk Tk; Ik]'*r
                    // ------------------------------------------------------------------------------------------------
                    s.noalias() = aRegularizationFactor * r.tail(ObjRank+RemainingColumns);

                    s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
                    s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
                    // ------------------------------------------------------------------------------------------------
                    gamma_previous = gamma;
                    gamma = s.squaredNorm();
                    beta = gamma/gamma_previous;
                    p = s + beta*p;
                    // ------------------------------------------------------------------------------------------------
                    iter++;
                }

                return iter;
            }

            /**
               \brief Accumulate nullspace basis.

               \note This implementation is counterintuitive.

               \verbatim
               For example, consider the first two null space bases
               Z1 = [-inv(R1)*T1; I] = [A,B;I], where A has as many columns as R2
               Z2 = [-inv(R2)*T2; I]

               So, essentially we have to perform Z1*Z2 = [-A*inv(R2)*T2 + B; Z2]

               We do this as follows:

               1. The "null_space" array before the multiplication has the form ("I" is not stored)

               [A, B; 0, 0]

               2. Initialize the lower left block

               [A, B; I, 0]

               3. Perform one triangular solve [A;I]*inv(R2) (note that inv(R2) is not negated)

               [A*inv(R2), B; inv(R2), 0]

               4. Multiply (note the sign)
               [B;0] = [B;0] - [A*inv(R2);inv(R2)]*T2
               \endverbatim

               Note that we don't exploit the fact that R2, and hence inv(R2), is upper
               triangular. In my tests, however this behaved faster than
               accumulate_nullspace_basis_1

               \todo the above observation is worth verifying again
            */
            void accumulate_nullspace_basis(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                dBlockType Rk(LOD,
                              FirstRowIndex, FirstColIndex,
                              ObjRank, ObjRank);

                dBlockType UpBlock(LOD,
                                   FirstRowIndex, FirstColIndex + ObjRank,
                                   ObjRank, RemainingColumns+1);

                dBlockType LeftBlock(null_space,
                                     0, FirstColIndex,
                                     FirstColIndex-nVarFixed+ObjRank, ObjRank);

                dBlockType TrailingBlock(null_space,
                                         0, FirstColIndex + ObjRank,
                                         FirstColIndex-nVarFixed+ObjRank, RemainingColumns+1);

                LeftBlock.block(FirstColIndex-nVarFixed,0,ObjRank,ObjRank).setIdentity();

                Rk.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(LeftBlock);

                if (ObjRank == 1)
                {
                    TrailingBlock.noalias() -= LeftBlock.col(0) * UpBlock.row(0);
                }
                else if (ObjRank > 2 && ObjRank <= 8)
                {
                    for (Index k=0; k<RemainingColumns+1; k++)
                    {
                        TrailingBlock.col(k).noalias() -= LeftBlock * UpBlock.col(k);
                    }
                }
                else if ((ObjRank == 2) || (ObjRank > 8))
                {
                    TrailingBlock.noalias() -= LeftBlock * UpBlock;
                }
            }

            /**
               \brief Accumulate nullspace basis.

               \verbatim
               For example, consider the first two null space bases
               Z1 = [-inv(R1)*T1; I] = [A,B;I], where A has as many columns as R2
               Z2 = [-inv(R2)*T2; I]

               So, essentially we have to perform Z1*Z2 = [-A*inv(R2)*T2 + B; Z2]

               We do this as follows:

               1. The "null_space" array before the multiplication has the form ("I" is not stored)
               [A, B; 0, 0]

               2. Initialize the lower right block
               [A, B; 0, T2]

               3. Perform one triangular solve C = -inv(R2)*T2
               [A, B; 0, C]

               4. Update block B
               B = B + A*C
               \endverbatim

               \note This function has not been tested properly (with simple bounds)
            */
            void accumulate_nullspace_basis_1(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
            {
                dBlockType Rk(LOD,
                              FirstRowIndex, FirstColIndex,
                              ObjRank, ObjRank);

                dBlockType Tk(LOD,
                              FirstRowIndex, FirstColIndex + ObjRank,
                              ObjRank, RemainingColumns+1);

                // A
                dBlockType LeftPrevious(null_space,
                                        0, FirstColIndex,
                                        FirstColIndex-nVarFixed, ObjRank);

                // B
                dBlockType RightPrevious(null_space,
                                         0, FirstColIndex + ObjRank,
                                         FirstColIndex-nVarFixed, RemainingColumns+1);

                // C
                dBlockType iRkTk(null_space,
                                 FirstColIndex-nVarFixed, FirstColIndex + ObjRank,
                                 ObjRank, RemainingColumns+1);

                // C = -inv(R2)*T2
                iRkTk = -Tk;
                Rk.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(iRkTk);

                // B = B + A*C
                if (ObjRank == 1)
                {
                    RightPrevious.noalias() += LeftPrevious.col(0) * iRkTk.row(0);
                }
                else if (ObjRank > 2 && ObjRank <= 8)
                {
                    for (Index k=0; k<RemainingColumns+1; k++)
                    {
                        RightPrevious.col(k).noalias() += LeftPrevious * iRkTk.col(k);
                    }
                }
                else if ((ObjRank == 2) || (ObjRank > 8))
                {
                    RightPrevious.noalias() += LeftPrevious * iRkTk;
                }
            }

            // ==================================================================
            // definition of scalars
            // ==================================================================

            /**
                \brief Number of variables (including fixed variables, see #nVarFixed)
            */
            Index nVar;

            /**
                \brief Number of fixed variables
            */
            Index nVarFixed;

            /**
                \brief Number of fixed variables that have been initialized

                \note It is used when variables to be fixed are defined one by one e.g.,
                #fixVariable(...).
            */
            Index nVarFixedInit;

            /**
                \brief Number of objectives
            */
            Index nObj;

            /**
               \brief Total number of constraints in the objectives of the problem (this does not
               include the constraints used to fix variables)
            */
            Index nCtr;

            /**
               \brief Rank of the hierarchy
            */
            Index TotalRank;

            /**
                \brief Regularization factor used at the current level

                \note Based on obj_info[:].regularization_factor
            */
            RealScalar aRegularizationFactor;

            // ==================================================================
            // definition of vectors of integers
            // ==================================================================

            /**
                \brief Stores the indexes of permuted columns
            */
            iVectorType column_permutations;

            /**
                \brief Indexes of fixed variables

                \attention This field is changed in #factorize().
            */
            iVectorType fixed_var_index;

            // ==================================================================
            // definition of vectors of doubles
            // ==================================================================

            /**
                \brief Workspace array of doubles

                \verbatim
                --------------------------------------
                factorize()
                --------------------------------------
                hhWorkspace     nVar+1
                ColNorms        nVar
                --------------------------------------
                Lagrange multipliers related
                --------------------------------------
                Lambda          maxObjDimSum + nVar
                tmp_RHS         maxObjDimSum

                the "+ nVar" in Lambda is necessary
                when there are fixed variables
                --------------------------------------
                Others
                --------------------------------------
                w               maxObjDimSum
                \endverbatim

                \note #LOD and #x are kept separately for convenience.
            */
            dVectorType dWorkspace;

            /**
                \brief A (primal) solution

                \note When TotalRank < nVar the problem is under-determined and x is a "basic solution"
                (see p. 106 of "Numerical Methods for Least Squares Problems" by Bjorck, 1996).
            */
            dVectorType x;

            /*
              \brief Householder scalars associated with the QR factorization for a (LexLSE) objective.

              \note computed during the factorization.
            */
            dVectorType hh_scalars;

            /*
              \brief Store the residuals when regularizing

              \note To be initialized with the right-hand-side vector
            */
            dVectorType residual_mu;


            // ==================================================================
            // definition of matrices
            // ==================================================================

            /**
                \brief Constraint matrix

                \note On input, LOD contains the matrix of constraints (the last column contains the RHS
                vector) of the lexicographic problem. On output, LOD contains the matrices L, Q and R
                from the LOD of the constraint matrix. The orthonormal matrices Q are
                stored in Householder factored form below the main diagonal of the leading triangular
                matrices of each objective, while the last column contains (LQ)^{-1}*RHS.

                \note The row-size of LOD is set to the the maximum number of envisioned constraints
                (excluding "fixing constraints"). This is usefull in the context of LexLSI where one
                instance of LexLSE is used to solve problems with different dimensoins.
            */
            dMatrixType LOD;

            /**
               \todo REMOVE (JUST TESTING)
            */
            dMatrixType PROBLEM_DATA;

            /**
                \brief nVar x (nVar+1) matrix containing the remaining null-space + right-hand-side

                \note Used for Tikhonov regularization or terminal objective of the form (x = 0)
            */
            dMatrixType null_space;

            /**
                \brief Used for temporary storage when doing Tikhonov regularization (option: A'*inv(A*A'+I)*b)
            */
            dMatrixType array;

            /**
               \brief Used in implementation of conjugate gradients (CG)
            */
            dMatrixType rqsp_work;

            /**
               \brief Used to store the solutions of all problems in the sequence (necessary for
               computing Lagrange multipliers when we regularize)
            */
            dMatrixType X_mu;

            dMatrixType X_mu_rhs;

            // ==================================================================
            // other
            // ==================================================================

            /**
                \brief Vector containing information about the constraints involved in each objective
            */
            std::vector<ObjectiveInfo> obj_info;

            /**
                \brief Type of fixed variables
            */
            std::vector<ConstraintActivationType> fixed_var_type;

            /**
                \brief Type of (active) constraints
            */
            std::vector<ConstraintActivationType> ctr_type;

            /**
                \brief Permutation matrix
            */
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;


            /// \brief Parameters of the solver.
            ParametersLexLSE parameters;

        }; // END class LexLSE

    } // END namespace internal

} // END namespace LexLS

#endif // LEXLSE
