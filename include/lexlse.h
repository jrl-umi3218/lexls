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

                Index dim = std::max(maxObjDimSum,nVar);
                dWorkspace.resize(2*dim + nVar + 1);

                column_permutations.resize(nVar);
            }

            /**
                \brief Factorize using Householder QR (with column pivoting) + classical implementation of Gauss elimination

                \note This is an unblocked implementation of both the QR decomposition (essentially,
                Eigen's ColPivHouseholderQR) and the Gauss elimination.

                \attention When the rank of the constraints involved in a given objective is equal to 1,
                many simplifications can be made. In the current implementation, this case is benefited
                from only when forming the Schur complement (during the Gauss elimination).
            */
            void factorize()
            {
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
                dMatrixType array;
                array.resize(nVar,nVar+1);
                array.setZero();

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
                dMatrixType array;
                array.resize(nVar,nVar+1);
                array.setZero();

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
               \brief Compute the solution that minimizes ||M.leftCols(nVar)*x - M.col(nVar)||^2.

               \note using the normal equations

               \note M is copied (because I modify it below)
            */
            void solveGeneralNorm(dMatrixType M)
            {
                dMatrixType array;
                array.resize(nVar,nVar+1);
                array.setZero();

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

                P.setIdentity();
                x.tail(nVar-nVarFixed).setZero(); // x.head(nVarFixed) has already been initialized in fixVariable(...)
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
