// Time-stamp: <2014-12-11 13:12:16 drdv>
#ifndef LEXLSE
#define LEXLSE

#include <utility.h>

namespace LexLS
{    
    /** 
        \brief Definition of a lexicographic least-squares problem with equality constraints

        \todo To describe how variables can be fixed.

        \todo Should I resize FixedVarType.resize(nVar) and FixedVarIndex.resize(nVar) only once (in
        active-set iterations)?

        \todo There are two tolerances that are fixed. see lexlse.findDescentDirection(...) and factorize()
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
            nVarFixedInit(0),
            LinearDependenceTolerance(1e-12),
            isFactorized(false), 
            regularizationType(REGULARIZATION_NONE),
            regularizationMaxIterCG(10),
            realSensitivityResidual(false),
            isSolved(false) {}
        
        /** 
            \param[in] nVar_   Number of variables (only number of elements in x, and not in the residuals w)
            \param[in] nObj_   Number of objectives
            \param[in] ObjDim_ Number of constraints involved in each objective
        */
        LexLSE(Index nVar_, Index nObj_, Index *ObjDim_):
            nVarFixed(0),
            nVarFixedInit(0),
            LinearDependenceTolerance(1e-12),
            regularizationType(REGULARIZATION_NONE),
            regularizationMaxIterCG(10),
            realSensitivityResidual(false),
            isFactorized(false), 
            isSolved(false)
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

            ObjInfo.resize(nObj);
            x.resize(nVar);
            P.resize(nVar);
            
            Index maxObjDimSum = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                maxObjDimSum += maxObjDim[ObjIndex];

            hh_scalars.resize(maxObjDimSum);
            CtrType.resize(maxObjDimSum, CONSTRAINT_TYPE_UNKNOWN);

            LQR.resize(maxObjDimSum,nVar+1); // store the RHS as well, thus the "+1"
            
            NullSpace.resize(nVar,nVar+1);
            array.resize(nVar,nVar+1);

            Index dim = std::max(maxObjDimSum,nVar);
            dWorkspace.resize(2*dim + nVar + 1);
            
            regularization.resize(nObj);
            ColPermutations.resize(nVar);

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
        */
        void factorize()
        { 
            /// @todo unused variable normTrailingRHS deleted
            RealScalar maxColNormValue, tau, PivotValue;
            Index RemainingRows, ObjRank, ObjDim, TotalRank, maxColNormIndex;

            Index RowIndex;                   // Current constraint
            Index FirstRowIndex;              // The same as ObjInfo[k].FirstRowIndex (already available)
            Index FirstColIndex;              // The same as ObjInfo[k].FirstColIndex (not yet computed)
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
                    coeff = FixedVarIndex.coeffRef(k);
                    ColPermutations.coeffRef(k) = coeff;
                    if (k != coeff) 
                        LQR.col(k).head(nCtr).swap(LQR.col(coeff).head(nCtr));

                    for (Index i=k+1; i<nVarFixed; i++)
                    {
                        if (FixedVarIndex.coeffRef(i) == k)
                        {
                            FixedVarIndex.coeffRef(i) = coeff; // FixedVarIndex is modified, but this is OK.
                            break;
                        }
                    }
                }                
                LQR.col(nVar).head(nCtr).noalias() -= LQR.block(0,0,nCtr,nVarFixed)*x.head(nVarFixed);
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
                    P.applyTranspositionOnTheRight(k, ColPermutations.coeff(k));
                
                isFactorized = true;
                                
                return; // early termination if possible
            }

            // ----------------------------------------------------
            dVectorBlockType    ColNorms(dWorkspace,    0, nVar  );
            dVectorBlockType hhWorkspace(dWorkspace, nVar, nVar+1);
            // ----------------------------------------------------

            // DEBUG: to remove
            DedicatedVariables.resize(nVar);            
            for (Index k=0; k<nVar; k++)
                DedicatedVariables.coeffRef(k) = k;
            
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++) // loop over all (explicitly defined) objectives
            {    
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex = ColIndex;
                ObjDim        = ObjInfo[ObjIndex].dim;

                for(Index k=ColIndex; k<nVar; k++) // initially compute the norms of the columns
                    ColNorms.coeffRef(k) = LQR.col(k).segment(FirstRowIndex,ObjDim).squaredNorm();

                // QR factorization of constraints involved in objective ObjIndex using the remaining variables 
                for (Index counter=0; counter<ObjDim; counter++) // loop over all constraints in a given objective
                {
                    RowIndex      = FirstRowIndex + counter; // current row to process
                    RemainingRows = ObjDim        - counter; // remaining rows in the current objective

                    maxColNormValue = ColNorms.tail(RemainingColumns).maxCoeff(&maxColNormIndex);
                    maxColNormIndex += ColIndex;
	
                    // the next two lines are for numerical stability 
                    // (I use them because this is what is done in ColPivHouseholderQR.h)
                    maxColNormValue = LQR.col(maxColNormIndex).segment(RowIndex,RemainingRows).squaredNorm();
                    ColNorms.coeffRef(maxColNormIndex) = maxColNormValue;
                    
                    // After we break, elimination is performed if there are objectives with lower priority and ObjInfo[ObjIndex].rank > 0
                    if (maxColNormValue < LinearDependenceTolerance)
                        break;

                    // --------------------------------------------------------------------------
                    // apply column permutation
                    // --------------------------------------------------------------------------
                    ColPermutations.coeffRef(ColIndex) = maxColNormIndex;
                    if(ColIndex != maxColNormIndex)
                    {
                        LQR.col(ColIndex).head(nCtr).swap(LQR.col(maxColNormIndex).head(nCtr));
                        std::swap(ColNorms.coeffRef(ColIndex), ColNorms.coeffRef(maxColNormIndex));

                        std::swap(DedicatedVariables.coeffRef(ColIndex), DedicatedVariables.coeffRef(maxColNormIndex)); // DEBUG

                        NullSpace.col(ColIndex).head(FirstColIndex).swap(NullSpace.col(maxColNormIndex).head(FirstColIndex)); // FOR REGULARIZATION
                    }

                    // --------------------------------------------------------------------------
                    // apply Householder transformations (on the RHS as well)
                    // --------------------------------------------------------------------------
                    // when RemainingRows = 1, since sqrt(maxColNormValue) >= LinearDependenceTolerance, the Householder matrix is the identity (tau = 0)
                    if (RemainingRows > 1) 
                    {
                        LQR.col(ColIndex).segment(RowIndex,RemainingRows).makeHouseholderInPlace(tau,PivotValue);
                        LQR.coeffRef(RowIndex,ColIndex) = PivotValue;
                        LQR.block(RowIndex,ColIndex+1,RemainingRows,RemainingColumns) // apply transformation on the RHS as well
                            .applyHouseholderOnTheLeft(LQR.col(ColIndex).segment(RowIndex+1,RemainingRows-1), 
                                                       tau, 
                                                       &hhWorkspace.coeffRef(0));
                        hh_scalars.coeffRef(FirstRowIndex+counter) = tau;
                    }
                    // --------------------------------------------------------------------------

                    ColIndex++;
                    RemainingColumns = nVar - ColIndex;

                    // terminate the QR factorization (after the elimination step below, the LQR factorization is terminated as well)
                    if (RemainingColumns == 0)
                        break; 

                    // update our table of squared norms of the columns 
                    // (note that above ColIndex is incremented and RemainingColumns is decreased)
                    if (RemainingRows > 0)
                        ColNorms.tail(RemainingColumns) -= LQR.row(RowIndex).segment(ColIndex,RemainingColumns).cwiseAbs2();

                } // END for (counter=0; counter<ObjDim; counter++)

                // Note that here ColIndex is the index of the next available variable

                ObjRank = ObjInfo[ObjIndex].rank = ColIndex - FirstColIndex; // store the rank

                // -----------------------------------------------------------------------
                // Regularization
                // -----------------------------------------------------------------------
                RealScalar conditioning_estimate;

                if (1) // constant damping factor
                {
                    damp_factor = regularization[ObjIndex];
                }
                else // variable damping factor (JUST TESTING, NOT READY YET)
                {
                    damp_factor = 0.0;
                    if (ObjRank > 0)
                    {                        
                        // -----------------------------------------------------------------------
                        // conditioninig estimation
                        // -----------------------------------------------------------------------
                        dVectorType rhs_tmp = LQR.col(nVar).segment(FirstRowIndex,ObjRank);

                        //print_eigen_matrix(rhs_tmp, "rhs");

                        conditioning_estimate = rhs_tmp.squaredNorm();
                        LQR.block(FirstRowIndex, FirstColIndex, ObjRank, ObjRank).
                            triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(rhs_tmp);

                        //print_eigen_matrix(rhs_tmp, "x");

                        conditioning_estimate /= rhs_tmp.squaredNorm();
                        // -----------------------------------------------------------------------
                    }
                }
                
                if (ObjRank > 0)
                {
                    switch(regularizationType){

                    case REGULARIZATION_TIKHONOV:

                        //printf("REGULARIZATION_TIKHONOV \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            if (FirstColIndex + ObjRank <= RemainingColumns)
                                regularize_tikhonov_2(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                            else
                                regularize_tikhonov_1(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }            
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;

                    case REGULARIZATION_TIKHONOV_CG:

                        //printf("REGULARIZATION_TIKHONOV_CG(%d) \n", regularizationMaxIterCG);

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_tikhonov_CG(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                            //regularize_tikhonov_CG_x0(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;
                            
                    case REGULARIZATION_R:

                        //printf("REGULARIZATION_TIKHONOV_R \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_R(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;

                    case REGULARIZATION_R_NO_Z:

                        //printf("REGULARIZATION_TIKHONOV_R_NO_Z \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_R_NO_Z(FirstRowIndex, FirstColIndex, ObjRank);
                        }
                        break;

                    case REGULARIZATION_RT_NO_Z:

                        //printf("REGULARIZATION_TIKHONOV_RT_NO_Z \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_RT_NO_Z(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        break;

                    case REGULARIZATION_RT_NO_Z_CG:

                        //printf("REGULARIZATION_RT_NO_Z_CG(%d) \n", regularizationMaxIterCG);

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_RT_NO_Z_CG(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        break;

                    case REGULARIZATION_TIKHONOV_1:

                        //printf("REGULARIZATION_TIKHONOV_1 \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_tikhonov_1(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;

                    case REGULARIZATION_TIKHONOV_2:

                        //printf("REGULARIZATION_TIKHONOV_2 \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_tikhonov_2(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;                            

                    case REGULARIZATION_TEST:

                        //printf("REGULARIZATION_TIKHONOV_RT_NO_Z_SMOOTH \n");

                        if ( !isEqual(damp_factor,0.0) ) 
                        {
                            regularize_RT_NO_Z_smooth(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        }
                        accumulate_nullspace_basis(FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);
                        break;

                    case REGULARIZATION_NONE:

                        //printf("REGULARIZATION_NONE \n");

                        // do nothing
                        break;
                    }
                }

                // -----------------------------------------------------------------------
                // Gauss transformation
                // -----------------------------------------------------------------------

                if (ObjIndex < nObj-1) // if there are objectives with lower priority
                {                    
                    if (ObjRank > 0) // if there are variables to eliminate 
                    {
                        FirstRowIndexNextObjective = FirstRowIndex + ObjDim; 
                        RemainingRows = nCtr-FirstRowIndexNextObjective; // VARIABLE REDEFINITION: remaining rows after the current objective
                        
                        // update the trailing block of the matrix LQR
                        // handle the RHS vector as well, hence the "+1" (recall that RemainingColumns = nVar-ColIndex)
                        // here, we cannot use directly LQR.bottomRightCorner(RemainingRows,RemainingColumns+1).noalias()

                        dBlockType LeftBlock(LQR,FirstRowIndexNextObjective, FirstColIndex, RemainingRows, ObjRank);
                        dBlockType UpBlock(LQR,FirstRowIndex,ColIndex,ObjRank,RemainingColumns+1);
                        dBlockType TrailingBlock(LQR,FirstRowIndexNextObjective,ColIndex,RemainingRows,RemainingColumns+1);
                        
                        LQR.block(FirstRowIndex,FirstColIndex,ObjRank,ObjRank)
                            .triangularView<Eigen::Upper>()
                            .solveInPlace<Eigen::OnTheRight>(LeftBlock);

                        if (ObjRank == 1)
                        {
                            // the .col(0) and .row(0) are important only for efficient computation
                            TrailingBlock.noalias() -= LeftBlock.col(0) * UpBlock.row(0);
                        }
                        else if (ObjRank >= 2 && ObjRank <= 8)
                        {
                            for (Index k=0; k<RemainingColumns+1; k++)
                                TrailingBlock.col(k).noalias() -= LeftBlock * UpBlock.col(k);
                        }                        
                        else if (ObjRank > 8)
                        {
                            TrailingBlock.noalias() -= LeftBlock * UpBlock;
                        }
                    }
                }

                // -----------------------------------------------------------------------

                if (RemainingColumns == 0) // ColIndex = nVar
                {
                    // Initialize some remaining fields of ObjInfo before we terminate (used later in ComputeLambda())
                    for (Index k=ObjIndex+1; k<nObj; k++) 
                        ObjInfo[k].FirstColIndex = ObjInfo[k-1].FirstColIndex + ObjInfo[k-1].rank; // of course, ObjInfo[>ObjIndex].rank = 0
                    
                    break; // terminate the LQR factorization
                }

            } // END for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)

            TotalRank = nVarFixed;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                TotalRank += ObjInfo[ObjIndex].rank;
          
            // form the permutation matrix
            for (Index k=0; k<TotalRank; k++)
                P.applyTranspositionOnTheRight(k, ColPermutations.coeff(k));

            isFactorized = true;

        } // END factorize()
        
        /** 
            \brief Computes the sensitivity of objective ObjIndex with respect to (small) variatoins
            of the constraints involved in the LexLSE problem

            \note Upon exit, the lagrange multipliers associated with objective ObjIndex can be accessed using
            dWorkspace.head(nVarFixed + nLambda)

            \todo Replace applyOnTheLeft(householderSequence(H,h)) with my implementation
        */
        bool ObjectiveSensitivity(Index ObjIndex, 
                                  Index &CtrIndex2Remove, Index &ObjIndex2Remove, 
                                  RealScalar tolWrongSignLambda, RealScalar tolCorrectSignLambda,
                                  dVectorType &residual)
        {
            RealScalar maxAbsValue = 0.0;
            bool tmp_bool, FoundBetterDescentDirection = false;
            Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
            Index &ColDim = FirstColIndex; // ColDim will be used to indicate column dimension (I use it for clarity)

            // Even though this computation can be improved, it is very cheap and it is convenient to perform here 
            Index nLambda = 0; // Number of constraints that may influence objective ObjIndex (excluding "variable fixing" constraints)
            Index nRank   = 0; // Total rank of objectives before objective ObjIndex
            for (Index k=0; k<ObjIndex; k++)
            {
                nLambda += ObjInfo[k].dim;
                nRank   += ObjInfo[k].rank;
            }
            nLambda += ObjInfo[ObjIndex].dim;
            
            // ---------------------------------------------------------------------------------
            dWorkspace.head(nVarFixed + nLambda + nRank + nVarFixed).setZero();
            dVectorBlockType LambdaFixed(dWorkspace,                   0, nVarFixed);
            dVectorBlockType      Lambda(dWorkspace,           nVarFixed, nLambda);
            dVectorBlockType         rhs(dWorkspace, nLambda + nVarFixed, nRank+nVarFixed);
            // ---------------------------------------------------------------------------------

            FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
            FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
            ObjDim        = ObjInfo[ObjIndex].dim;
            ObjRank       = ObjInfo[ObjIndex].rank;

            // Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed

            if (realSensitivityResidual) // use the real residual (not recommended)
            {
                Lambda.segment(FirstRowIndex, ObjDim) = residual.head(ObjDim);
            }
            else // compute the residual from the factorization
            {
                // copy only what is needed to compute the residual w = A*x-b (i.e., -y_hat)
                Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = \
                    -LQR.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);
                
                // compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
                Lambda.segment(FirstRowIndex, ObjDim)
                    .applyOnTheLeft(householderSequence(LQR.block(FirstRowIndex, 
                                                                  FirstColIndex, 
                                                                  ObjDim, 
                                                                  ObjRank),
                                                        hh_scalars.segment(FirstRowIndex,ObjDim))); 
            }
            
            // check for wrong sign of the Lagrange multipliers
            FoundBetterDescentDirection = findDescentDirection(FirstRowIndex,
                                                               ObjDim,
                                                               maxAbsValue,
                                                               CtrIndex2Remove,
                                                               Lambda,
                                                               tolWrongSignLambda,
                                                               tolCorrectSignLambda);
            
            if (FoundBetterDescentDirection)
                ObjIndex2Remove = ObjIndex;  
            
            if (ObjIndex>0) // the first objective has only Lagrange multipliers equal to the optimal residual   
            {
                // e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}
                rhs.head(ColDim).noalias() = -LQR.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                    Lambda.segment(FirstRowIndex, ObjDim);
                
                for (Index k=ObjIndex-1; k>=0; k--)
                {
                    FirstRowIndex = ObjInfo[k].FirstRowIndex;
                    FirstColIndex = ObjInfo[k].FirstColIndex;
                    ObjDim        = ObjInfo[k].dim;
                    ObjRank       = ObjInfo[k].rank;

                    // Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed
                    
                    Lambda.segment(FirstRowIndex, ObjRank) = rhs.segment(FirstColIndex, ObjRank);
                    
                    // apply Q_k' on the left
                    Lambda.segment(FirstRowIndex, ObjDim)
                        .applyOnTheLeft(householderSequence(LQR.block(FirstRowIndex, 
                                                                      FirstColIndex, 
                                                                      ObjDim, 
                                                                      ObjRank),
                                                            hh_scalars.segment(FirstRowIndex,ObjDim)));
                    
                    rhs.head(ColDim).noalias() -= LQR.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                        Lambda.segment(FirstRowIndex, ObjDim);
                    
                    // check for wrong sign of the Lagrange multipliers
                    tmp_bool = findDescentDirection(FirstRowIndex,
                                                    ObjDim,
                                                    maxAbsValue,
                                                    CtrIndex2Remove,
                                                    Lambda,
                                                    tolWrongSignLambda,
                                                    tolCorrectSignLambda);
                    
                    if (tmp_bool)
                        ObjIndex2Remove = k;
                    FoundBetterDescentDirection = (FoundBetterDescentDirection || tmp_bool);
                } // END for (Index k=ObjIndex-1; k>=0; k--)
            }
            
            if (nVarFixed>0) // Handle fixed variables (if any)
            {
                LambdaFixed = -LQR.block(0, 0, nLambda, nVarFixed).transpose() * Lambda;

                // check for wrong sign of the Lagrange multipliers
                tmp_bool = findDescentDirection(-1,
                                                nVarFixed,
                                                maxAbsValue,
                                                CtrIndex2Remove,
                                                LambdaFixed,
                                                tolWrongSignLambda,
                                                tolCorrectSignLambda);
                                
                // -1(-st) objective implies the fixed variables
                //   :if there are no fixed variables we will not be here
                //   :if there are fixed variables ObjIndex2Remove + ObjOffset is used in LexLSI
                if (tmp_bool)
                    ObjIndex2Remove = -1;
                FoundBetterDescentDirection = (FoundBetterDescentDirection || tmp_bool);
            }
/*
            if (FoundBetterDescentDirection)
            {
                printf("-------------------------------------------- gamma = %+e\n", maxAbsValue);
            }
*/
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
                nLambda += ObjInfo[k].dim;
                nRank   += ObjInfo[k].rank;
            }
            nLambda += ObjInfo[ObjIndex].dim;
            
            // ---------------------------------------------------------------------------------
            dWorkspace.head(nVarFixed + nLambda + nRank + nVarFixed).setZero();
            dVectorBlockType LambdaFixed(dWorkspace,                   0, nVarFixed);
            dVectorBlockType      Lambda(dWorkspace,           nVarFixed, nLambda);
            dVectorBlockType         rhs(dWorkspace, nLambda + nVarFixed, nRank+nVarFixed);
            // ---------------------------------------------------------------------------------

            FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
            FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
            ObjDim        = ObjInfo[ObjIndex].dim;
            ObjRank       = ObjInfo[ObjIndex].rank;

            // Lambda.segment(FirstRowIndex, ObjRank).setZero(); assumed
            
            // copy only what is needed to compute the residual w = A*x-b (i.e., -y_hat)
            Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) =     \
                -LQR.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);
            
            // compute the optimal residual associated with objective ObjIndex (apply Q_{ObjIndex} on the left)
            Lambda.segment(FirstRowIndex, ObjDim)
                .applyOnTheLeft(householderSequence(LQR.block(FirstRowIndex, 
                                                              FirstColIndex, 
                                                              ObjDim, 
                                                              ObjRank),
                                                    hh_scalars.segment(FirstRowIndex,ObjDim))); 
                        
            if (ObjIndex>0) // the first objective has only Lagrange multipliers equal to the optimal residual   
            {
                // e.g., for the fourth objective, here we perform [L41, L42, L43]' * {optimal residual from above}
                rhs.head(ColDim).noalias() = -LQR.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                    Lambda.segment(FirstRowIndex, ObjDim);
                
                for (Index k=ObjIndex-1; k>=0; k--)
                {
                    FirstRowIndex = ObjInfo[k].FirstRowIndex;
                    FirstColIndex = ObjInfo[k].FirstColIndex;
                    ObjDim        = ObjInfo[k].dim;
                    ObjRank       = ObjInfo[k].rank;
                    
                    // Lambda.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank).setZero(); assumed
                    
                    Lambda.segment(FirstRowIndex, ObjRank) = rhs.segment(FirstColIndex, ObjRank);
                    
                    // apply Q_k' on the left
                    Lambda.segment(FirstRowIndex, ObjDim)
                        .applyOnTheLeft(householderSequence(LQR.block(FirstRowIndex, 
                                                                      FirstColIndex, 
                                                                      ObjDim, 
                                                                      ObjRank),
                                                            hh_scalars.segment(FirstRowIndex,ObjDim)));
                    
                    rhs.head(ColDim).noalias() -= LQR.block(FirstRowIndex, 0, ObjDim, ColDim).transpose() * \
                        Lambda.segment(FirstRowIndex, ObjDim);

                } // END for (Index k=ObjIndex-1; k>=0; k--)
            }
            
            if (nVarFixed>0) // Handle fixed variables (if any)
                LambdaFixed = -LQR.block(0, 0, nLambda, nVarFixed).transpose() * Lambda;
            
            rhs.setZero(); // for convenience (easier to analyze the Lagrange multipliers by hand) 

        } // END ObjectiveSensitivity(...)
        
        /**
           \brief Given a vector of Lagrange multipliers, determine the largest (in absolute value)
           multiplier with a wrong sign.

           \param[in]     FirstRowIndex  Index of first element of the current objective (if FirstRowIndex < 0 handle fixed variables)
           \param[in]     ObjDim         Number of constraints.
           \param[in,out] maxAbsValue    Largest (in absolute value) multiplier with wrong sign
           \param[in,out] CtrIndex       Index of largest (in absolute value) multiplier with wrong sign.
           \param[in]     lambda         Vector of lagrange multipliers.
           \param[in]     tolWrongSignLambda   Absolute value of Lagrange multiplier to be considered with "wrong" sign.
           \param[in]     tolCorrectSignLambda Absolute value of Lagrange multiplier to be considered with "correct" sign.

           \note Lagrange multipliers in the interval (-tolWrongSignLambda --- 0 --- tolCorrectSignLambda) are considered equal to zero.

           \note Using tolCorrectSignLambda = 0 is not a good idea in general. This might lead to
           situations where in order to decrease the norm of the residual of a higher-level task
           with e.g., 1e-12, the solver might worsen the nor of the residual of lower-level taks
           with a lot more (e.g., 10).

           \return true if there are multipliers with a wrong sign whose absolute value is larger
           than the largest multiplier with a wrong sign from previous groups of Lagrange
           multipliers.
        */
        bool findDescentDirection(Index FirstRowIndex,
                                  Index ObjDim,
                                  RealScalar &maxAbsValue,   // modified
                                  Index &CtrIndex,           // modified
                                  const dVectorType& lambda,
                                  RealScalar tolWrongSignLambda, 
                                  RealScalar tolCorrectSignLambda)
        {
            bool FoundBetterDescentDirection = false;
            RealScalar aLambda;
            Index ind;
            ConstraintType *aCtrType;

            for (Index k=0; k<ObjDim; k++)
            {
                if (FirstRowIndex < 0) // handle fixed variables
                {
                    ind      = k;
                    aCtrType = &FixedVarType[ind];
                }
                else                   // handle general constraints
                {           
                    ind      = FirstRowIndex+k;
                    aCtrType = &CtrType[ind];
                }

                if (*aCtrType != EQUALITY_CONSTRAINT && *aCtrType != CORRECT_SIGN_OF_LAMBDA)
                {
                    aLambda = lambda.coeff(ind);

                    if (*aCtrType == LOWER_BOUND)
                        aLambda = -aLambda;
                            
                    // FIXME: to have as user input
                    if (aLambda > tolCorrectSignLambda) // is this reasonable?
                    {
                        *aCtrType = CORRECT_SIGN_OF_LAMBDA;
                    }
                    else if (aLambda < -tolWrongSignLambda)
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
            for(Index k=nObj-1; k>=0; k--) 
            {
                ObjRank = ObjInfo[k].rank;
                if (ObjRank > 0)
                {
                    dVectorBlockType x_k(x, ObjInfo[k].FirstColIndex, ObjRank);

                    x_k = LQR.col(nVar).segment(ObjInfo[k].FirstRowIndex, ObjRank);
                    
                    if (AccumulatedRanks > 0) // Do not enter here the first time ObjRank != 0
                    {
                        x_k.noalias() -= LQR.block(ObjInfo[k].FirstRowIndex,
                                                   ObjInfo[k+1].FirstColIndex,
                                                   ObjRank,
                                                   AccumulatedRanks) * x.segment(ObjInfo[k+1].FirstColIndex, AccumulatedRanks);
                    }
                    
                    LQR.block(ObjInfo[k].FirstRowIndex, ObjInfo[k].FirstColIndex, ObjRank, ObjRank)
                        .triangularView<Eigen::Upper>().solveInPlace(x_k);
                    
                    AccumulatedRanks += ObjRank;
                }
            }
            
            // Apply permutation
            x = P*x;

            // Problem solved
            isSolved = true;
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
                nVarRank += ObjInfo[ObjIndex].rank;           

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
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;
                
                RT.block(counter, counter, ObjRank, col_dim)
                    .triangularView<Eigen::Upper>() = LQR.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim);
                
                rhs.segment(counter,ObjRank) = LQR.col(nVar).segment(FirstRowIndex,ObjRank);

                counter += ObjRank;
                col_dim -= ObjRank;
            }

            // -------------------------------------------------------------------------
            // zero-out the redundant part (by applying Givens rotations on the right)
            // -------------------------------------------------------------------------
            GivensRotationSequence gs(nVarFree*nVarRank);
            for (Index i=0; i<nVarFree; i++)
            {
                for (Index j=nVarRank-1; j>=0; j--)
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
            for (Index i=gs.size()-1; i>=0; i--)
                rhs.applyOnTheLeft(gs.get_i(i), gs.get_j(i), gs.get(i));

            // -------------------------------------------------------------------------
            // Apply permutation
            // -------------------------------------------------------------------------
            x.tail(nVarRank + nVarFree) = rhs; // x.head(nVarFixed) contain the fixed variables
            x = P*x;

            // Problem solved
            isSolved = true;
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
                nVarRank += ObjInfo[ObjIndex].rank;           

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
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;
              
                RT.block(counter, counter, ObjRank, col_dim+1)
                    .triangularView<Eigen::Upper>() = LQR.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim+1);

                counter += ObjRank;
                col_dim -= ObjRank;
            }

            // no need to negate
            R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(T);

            D.triangularView<Eigen::Lower>() = T.leftCols(nVarFree).transpose()*T.leftCols(nVarFree);
            for (Index i=0; i<nVarFree; i++)
                D.coeffRef(i,i) += 1.0;

            d = T.leftCols(nVarFree).transpose() * T.col(nVarFree);

            Eigen::LLT<MatrixType> chol(D);
            x.tail(nVarFree) = chol.solve(d);

            counter = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;

                x.segment(nVarFixed+counter,ObjRank) = LQR.col(nVar).segment(FirstRowIndex,ObjRank) - 
                    LQR.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);
                                            
                counter += ObjRank;
            }
            R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(x.segment(nVarFixed,nVarRank));

            // -------------------------------------------------------------------------
            // Apply permutation
            // -------------------------------------------------------------------------
            x = P*x;

            // Problem solved
            isSolved = true;
        }        

        /**
           \brief Compute the least-norm solution.

           \note using the normal equations (and reusing the basis constructed for Tikhonov
           regularization). Hence in order to use this, one has to set regularizationType =
           REGULARIZATION_TIKHONOV and regularization = zeros(1,nObj);
        */
        void solveLeastNorm_3()
        {
            Index FirstRowIndex, FirstColIndex, ObjRank, counter = 0;
            Index nVarRank = 0; // number of variables determined by rank([R,T])

            // -------------------------------------------------------------------------
            // determine dimensions
            // -------------------------------------------------------------------------
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                nVarRank += ObjInfo[ObjIndex].rank;           

            Index nVarFree = nVar - (nVarRank + nVarFixed);

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType iR(NullSpace, nVarFixed,          nVarFixed, nVarRank,   nVarRank); // inv(R)
            dBlockType  T(NullSpace, nVarFixed, nVarFixed+nVarRank, nVarRank, nVarFree+1); // inv(R)*T
            dBlockType  D(    array,         0,                  0, nVarFree,   nVarFree);

            dBlockType2Vector d(array, 0, nVar, nVarFree, 1);

            // -------------------------------------------------------------------------

            D.triangularView<Eigen::Lower>() = T.leftCols(nVarFree).transpose()*T.leftCols(nVarFree);
            for (Index i=0; i<nVarFree; i++)
                D.coeffRef(i,i) += 1.0;

            d = T.leftCols(nVarFree).transpose() * T.col(nVarFree);

            Eigen::LLT<MatrixType> chol(D);
            x.tail(nVarFree) = chol.solve(d);

            counter = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;

                x.segment(nVarFixed+counter,ObjRank) = LQR.col(nVar).segment(FirstRowIndex,ObjRank) - 
                    LQR.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);
                                            
                counter += ObjRank;
            }
            x.segment(nVarFixed,nVarRank) = iR.triangularView<Eigen::Upper>() * x.segment(nVarFixed,nVarRank);

            // -------------------------------------------------------------------------
            // Apply permutation
            // -------------------------------------------------------------------------
            x = P*x;

            // Problem solved
            isSolved = true;
        }


        /**
           \brief Compute the solution that minimizes ||M.leftCols(nVar)*x - M.col(nVar)||^2.

           \note using the normal equations

           \note M is copied (because I modify it below) 
        */
        void solveGeneralNorm(MatrixType M)
        {
            Index FirstRowIndex, FirstColIndex, ObjRank, counter = 0;
            Index nVarRank = 0; // number of variables determined by rank([R,T])

            M = M*P; // permute columns

            // -------------------------------------------------------------------------
            // determine dimensions
            // -------------------------------------------------------------------------
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                nVarRank += ObjInfo[ObjIndex].rank;           

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
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;
              
                RT.block(counter, counter, ObjRank, col_dim+1)
                    .triangularView<Eigen::Upper>() = LQR.block(FirstRowIndex,FirstColIndex,ObjRank,col_dim+1);

                counter += ObjRank;
                col_dim -= ObjRank;
            }

            // -------------------------------------------------------------------------

            R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(LB);
            TB.noalias() -= LB * T;

            D.triangularView<Eigen::Lower>() = TB.leftCols(nVarFree).transpose()*TB.leftCols(nVarFree);
            d = TB.leftCols(nVarFree).transpose() * TB.col(nVarFree);

            Eigen::LLT<MatrixType> chol(D);
            x.tail(nVarFree) = chol.solve(d);

            counter = 0;
            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjRank       = ObjInfo[ObjIndex].rank;

                x.segment(nVarFixed+counter,ObjRank) = LQR.col(nVar).segment(FirstRowIndex,ObjRank) - 
                    LQR.block(FirstRowIndex,nVarRank+nVarFixed,ObjRank,nVarFree)*x.tail(nVarFree);
                                            
                counter += ObjRank;
            }
            R.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheLeft>(x.segment(nVarFixed,nVarRank));

            // -------------------------------------------------------------------------
            // Apply permutation
            // -------------------------------------------------------------------------
            x = P*x;

            // Problem solved
            isSolved = true;
        }        

        // =================================================================================================
        // set & get
        // =================================================================================================

        /**
           \brief Declare a variable as fixed

           \param[in] VarIndex Index of variable: x(VarIndex) will be fixed
           \param[in] VarValue Fix x(VarIndex) = VarValue
           \param[in] type     LOWER_BOUND if x(VarIndex) is fixed at its lower bound, UPPER_BOUND if
                               x(VarIndex) is fixed at its upper bound in LexLSI

           \note When we don't care about the type (i.e., we will not use the associated Lagrange
           multipliers), UPPER_BOUND is (arbitrary) chosen as a default value for type.

        */
        void fixVariable(Index VarIndex, RealScalar VarValue, ConstraintType type = UPPER_BOUND)
        {
            FixedVarIndex(nVarFixedInit) = VarIndex;
            x(nVarFixedInit)             = VarValue;
            FixedVarType[nVarFixedInit]  = type;
            
            nVarFixedInit++;
        }

        /** 
            \brief Declare some of the variables as fixed

            \param[in] nVarFixed_ Number of fixed variables
            \param[in] VarIndex   Indexes of variables to fix
            \param[in] VarValue   Values of the fixed variables
            \param[in] type       can be LOWER_BOUND or UPPER_BOUND
        */
        void fixVariables(Index nVarFixed_, Index *VarIndex, RealScalar *VarValue, ConstraintType *type)
        {
            if (nVarFixed_ > nVar)
                throw Exception("Cannot fix more than nVar variables");
            else
                nVarFixed = nVarFixed_;

            FixedVarIndex.resize(nVarFixed);
            FixedVarType.resize(nVarFixed);
                
            // copy data
            FixedVarIndex     = Eigen::Map<iVectorType>(VarIndex, nVarFixed);
            x.head(nVarFixed) = Eigen::Map<dVectorType>(VarValue, nVarFixed); // x will be permuted later
            FixedVarType.assign(type, type+nVarFixed);
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

                ObjInfo[ObjIndex].dim = ObjDim_[ObjIndex];
                if (ObjIndex > 0)
                    ObjInfo[ObjIndex].FirstRowIndex = ObjInfo[ObjIndex-1].FirstRowIndex + ObjInfo[ObjIndex-1].dim;
            }
            
            initialize();
        }

        /** 
            \brief Set the tolerances

            \param[in] LinearDependenceTolerance_ #LinearDependenceTolerance
        */      
        void setTolerance(RealScalar LinearDependenceTolerance_)
        {
            LinearDependenceTolerance = LinearDependenceTolerance_;
        }

        /** 
            \brief Set number of fixed variables

            \param[in] nVarFixed_ Number of fixed variables
        */        
        void setFixedVariablesCount(Index nVarFixed_)
        {
            if (nVarFixed_ > nVar)
                throw Exception("Cannot fix more than nVar variables");
            else
                nVarFixed = nVarFixed_;
            
            FixedVarIndex.resize(nVarFixed);
            FixedVarType.resize(nVarFixed);
        }

        /** 
            \brief Set (a non-negative) regularization factor for objective ObjIndex

            \note Note that fixed variables are not counted as an objective
        */
        void setRegularization(Index ObjIndex, RealScalar RegularizationFactor)
        {
            // @todo: check whether ObjIndex and RegularizationFactor make sense. 

            regularization(ObjIndex) = RegularizationFactor;
        }

        void setRegularizationType(RegularizationType regularizationType_)
        {
            regularizationType = regularizationType_;
        }

        void setRegularizationMaxIterCG(Index regularizationMaxIterCG_)
        {
            regularizationMaxIterCG = regularizationMaxIterCG_;
        }
        
        void setRealSensitivityResidual(bool realSensitivityResidual_)
        {
            realSensitivityResidual = realSensitivityResidual_;
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
            return FixedVarIndex;
        }

        /** 
            \brief set a random problem [A,RHS]
        */
        void setProblem(const MatrixType& data)
        {
            LQR = data;
        }

        /** 
            \brief Set data of objective ObjIndex [A,RHS] - the data is copied

            \param[in] ObjIndex Index of objective
            \param[in] data     data (including LHS & RHS)
        */
        void setData(Index ObjIndex, const MatrixType& data)
        {
            if (ObjIndex >= nObj)                
                throw Exception("ObjIndex >= nObj");
            
            LQR.block(ObjInfo[ObjIndex].FirstRowIndex,0,ObjInfo[ObjIndex].dim,nVar+1) = data;
        }

        /** 
            \brief Set one constraint

            \param[in] CtrIndex Index of row in LQR (regardless of objective)
            \param[in] row      Constraint vector (LHS)
            \param[in] rhs      RHS vector
        */
        void setCtr(Index CtrIndex, const dRowVectorType& row, RealScalar rhs)
        {
            LQR.row(CtrIndex).head(nVar) = row;
            LQR.coeffRef(CtrIndex,nVar)  = rhs;
        }
        
        /** 
            \brief Set the type of a constraint with index CtrIndex in LexLSE objective ObjIndex
        */      
        void setCtrType(Index ObjIndex, Index CtrIndex, ConstraintType type)
        {
            Index FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
            CtrType[FirstRowIndex + CtrIndex] = type;
        }
        
        /**
           \brief Form the residuals (A*x-RHS) through the LQR factorizetion. The residual of the
           fixed variables is always zero (and is not included).

           \note This function could compute an incorect residual depending on #LinearDependenceTolerance.
        */
        dVectorType& getResidual()
        {
            Index FirstRowIndex, FirstColIndex, ObjDim, ObjRank;
            dVectorBlockType w(dWorkspace,0,nCtr);

            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
            {
                FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
                FirstColIndex = ObjInfo[ObjIndex].FirstColIndex;
                ObjDim        = ObjInfo[ObjIndex].dim;
                ObjRank       = ObjInfo[ObjIndex].rank;
                
                w.segment(FirstRowIndex, ObjRank).setZero(); // Zero-out first ObjRank elements
                w.segment(FirstRowIndex+ObjRank, ObjDim-ObjRank) = -LQR.col(nVar).segment(FirstRowIndex+ObjRank, ObjDim-ObjRank);
                
                w.segment(FirstRowIndex, ObjDim)
                    .applyOnTheLeft(householderSequence(LQR.block(FirstRowIndex, 
                                                                  FirstColIndex, 
                                                                  ObjDim, 
                                                                  ObjRank),
                                                        hh_scalars.segment(FirstRowIndex,ObjDim)));
            }
                        
            return dWorkspace; // return directly a reference to the working array (use dWorkspace.head(nCtr) outside)
        }

        /** 
            \brief Return the type of constraint CtrIndex in objective ObjIndex
        */        
        ConstraintType getCtrType(Index ObjIndex, Index CtrIndex)
        {
            Index FirstRowIndex = ObjInfo[ObjIndex].FirstRowIndex;
            return CtrType[FirstRowIndex+CtrIndex];
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
            return ObjInfo[ObjIndex].dim;
        }

        /** 
            \brief Return the rank of the equations in objective ObjIndex
        */        
        Index getRank(Index ObjIndex) const
        {
            return ObjInfo[ObjIndex].rank;
        }

        /** 
            \brief Return the constraint matrix
        */
        MatrixType& getLQR()
        {
            return LQR;
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
            isFactorized  = false;
            isSolved      = false;

            for (Index ObjIndex=0; ObjIndex<nObj; ObjIndex++)
                ObjInfo[ObjIndex].rank = 0; 

            hh_scalars.setZero();

            NullSpace.setIdentity(); //todo: I could set it to 0
            array.setZero();

            P.setIdentity();
            x.tail(nVar-nVarFixed).setZero(); // x.head(nVarFixed) has already been initialized in fixVariable(...)

            regularization.setZero(); // by default there is no regularization
        }

        /** 
            \brief Tikhonov regularization (using the normal equations inv(A'*A+I)*A'*b)

            \note fast when column-dimension is small

            todo: maybe use "array" instaed of NullSpace for temporary storage (this would impact accumulate_nullspace_basis)
        */        
        void regularize_tikhonov_1(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(      LQR, FirstRowIndex,         FirstColIndex,                  ObjRank,                  ObjRank);
            dBlockType Tk(      LQR, FirstRowIndex, FirstColIndex+ObjRank,                  ObjRank,         RemainingColumns);
            dBlockType up(NullSpace,             0,         FirstColIndex,            FirstColIndex, RemainingColumns+ObjRank);
            dBlockType  D(    array,             0,                     0, RemainingColumns+ObjRank, RemainingColumns+ObjRank);

            dBlockType2Vector d(array, 0, nVar, RemainingColumns+ObjRank, 1);

            // -------------------------------------------------------------------------

            D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
            D.block(ObjRank,ObjRank,RemainingColumns,RemainingColumns).triangularView<Eigen::Lower>() = Tk.transpose()*Tk;
            D.block(ObjRank,0,RemainingColumns,ObjRank).noalias() = Tk.transpose()*Rk.triangularView<Eigen::Upper>();
            D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;

            for (Index i=0; i<RemainingColumns+ObjRank; i++)
                D.coeffRef(i,i) += mu;

            // ==============================================================================================
            d.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);
            d.tail(RemainingColumns).noalias() = Tk.transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            d.noalias() += mu * up.transpose() * NullSpace.col(nVar).head(FirstColIndex);

            // ==============================================================================================

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>() * d.head(ObjRank);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk * d.tail(RemainingColumns);
        }

        /** 
            \brief Tikhonov regularization (option: A'*inv(A*A'+I)*b)

            \note fast when row-dimension is small
        */        
        void regularize_tikhonov_2(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(      LQR, FirstRowIndex,         FirstColIndex,               ObjRank,                    ObjRank);
            dBlockType Tk(      LQR, FirstRowIndex, FirstColIndex+ObjRank,               ObjRank,           RemainingColumns);
            dBlockType up(NullSpace,             0,         FirstColIndex,         FirstColIndex, RemainingColumns + ObjRank);
            dBlockType  D(    array,             0,                     0, FirstColIndex+ObjRank,      FirstColIndex+ObjRank);

            dBlockType2Vector d(array, 0, nVar, FirstColIndex+ObjRank, 1);
            // -------------------------------------------------------------------------

            // forming D is very expensive
            D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>()  = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
            D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>() += Tk*Tk.transpose();

            D.block(ObjRank,ObjRank,FirstColIndex,FirstColIndex).triangularView<Eigen::Lower>() = mu*up*up.transpose();

            D.block(ObjRank, 0, FirstColIndex, ObjRank).noalias()  = 
                damp_factor*up.leftCols(ObjRank)*Rk.triangularView<Eigen::Upper>().transpose();
            
            D.block(ObjRank, 0, FirstColIndex, ObjRank).noalias() += 
                damp_factor*up.rightCols(RemainingColumns)*Tk.transpose();
            
            for (Index i=0; i<FirstColIndex+ObjRank; i++)
                D.coeffRef(i,i) += mu;

            d.head(ObjRank) = LQR.col(nVar).segment(FirstRowIndex,ObjRank);
            d.segment(ObjRank,FirstColIndex) = damp_factor*NullSpace.col(nVar).head(FirstColIndex);

            // -------------------------------------------------------------------------

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);

            // -------------------------------------------------------------------------

            for (Index i=0; i<FirstColIndex+ObjRank; i++)
                D.coeffRef(i,i) -= mu;

            d = D.selfadjointView<Eigen::Lower>() * d;
            LQR.col(nVar).segment(FirstRowIndex,ObjRank) = d.head(ObjRank);
        }

        /** 
            \brief Tikhonov regularization (basic variables)
        */        
        void regularize_R(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType  Rk(      LQR, FirstRowIndex, FirstColIndex,       ObjRank, ObjRank);
            dBlockType  up(NullSpace,             0, FirstColIndex, FirstColIndex, ObjRank);
            dBlockType   D(    array,             0,             0,       ObjRank, ObjRank);

            dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

            // -------------------------------------------------------------------------

            D.triangularView<Eigen::Lower>()  = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
            D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;
            for (Index i=0; i<ObjRank; i++)
                D.coeffRef(i,i) += mu;

            d.noalias()  = up.transpose() * NullSpace.col(nVar).head(FirstColIndex);
            d *= mu;
            d.noalias() += Rk.triangularView<Eigen::Upper>().transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            // -------------------------------------------------------------------------

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = Rk.triangularView<Eigen::Upper>() * d;
        }

        /** 
            \brief Tikhonov regularization (basic variables no Z)
        */        
        void regularize_R_NO_Z(Index FirstRowIndex, Index FirstColIndex, Index ObjRank)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk( LQR, FirstRowIndex, FirstColIndex, ObjRank, ObjRank);
            dBlockType  D(array,            0,             0, ObjRank, ObjRank);

            dBlockType2Vector d(array, 0, nVar, ObjRank, 1);
            // -------------------------------------------------------------------------

            D.triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
            for (Index i=0; i<ObjRank; i++)
                D.coeffRef(i,i) += mu;

            d.noalias() = Rk.triangularView<Eigen::Upper>().transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            // -------------------------------------------------------------------------

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = Rk.triangularView<Eigen::Upper>() * d;
        }

        /** 
            \brief RT regularization (no Z): [R,T;I]*x = [b;0]
        */        
        void regularize_RT_NO_Z(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(  LQR, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
            dBlockType Tk(  LQR, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);
            dBlockType  D(array,             0,                     0, ObjRank,          ObjRank);
            
            dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

            // -------------------------------------------------------------------------

            D.triangularView<Eigen::Lower>() = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
            D.triangularView<Eigen::Lower>() += Tk*Tk.transpose();

            for (Index i=0; i<ObjRank; i++)
                D.coeffRef(i,i) += mu;

            d = LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);

            // -------------------------------------------------------------------------

            for (Index i=0; i<ObjRank; i++)
                D.coeffRef(i,i) -= mu;
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() = D.selfadjointView<Eigen::Lower>() * d;
        }

        void regularize_RT_NO_Z_smooth(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(      LQR, FirstRowIndex,         FirstColIndex,                  ObjRank,                  ObjRank);
            dBlockType Tk(      LQR, FirstRowIndex, FirstColIndex+ObjRank,                  ObjRank,         RemainingColumns);
            dBlockType up(NullSpace,             0,         FirstColIndex,            FirstColIndex, RemainingColumns+ObjRank);
            dBlockType  D(    array,             0,                     0, RemainingColumns+ObjRank, RemainingColumns+ObjRank);

            dBlockType2Vector d(array, 0, nVar, RemainingColumns+ObjRank, 1);

            MatrixType V(FirstColIndex, RemainingColumns+ObjRank);
            V.setZero();
            V.block(0,0,FirstColIndex,FirstColIndex).setIdentity();

            // -------------------------------------------------------------------------

            D.setZero();
            D.block(0,0,ObjRank,ObjRank).triangularView<Eigen::Lower>() = (Rk.transpose()*Rk.triangularView<Eigen::Upper>()).eval();
            D.block(ObjRank,ObjRank,RemainingColumns,RemainingColumns).triangularView<Eigen::Lower>() = Tk.transpose()*Tk;
            D.block(ObjRank,0,RemainingColumns,ObjRank).noalias() = Tk.transpose()*Rk.triangularView<Eigen::Upper>();
            //D.triangularView<Eigen::Lower>() += mu*up.transpose()*up;
            D.triangularView<Eigen::Lower>() += mu*MatrixType::Identity(RemainingColumns+ObjRank,RemainingColumns+ObjRank);
            //D.triangularView<Eigen::Lower>() += mu*V.transpose()*V;

            for (Index i=0; i<RemainingColumns+ObjRank; i++)
                D.coeffRef(i,i) += mu;
            
            // ==============================================================================================
            d.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);
            d.tail(RemainingColumns).noalias() = Tk.transpose()*LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            //d.noalias() += mu * up.transpose() * NullSpace.col(nVar).head(FirstColIndex);

            //d.noalias() += mu * V.transpose() * NullSpace.col(nVar).head(FirstColIndex);
            
            // ==============================================================================================

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>() * d.head(ObjRank);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk * d.tail(RemainingColumns);
        }

        /** 
            \brief Tikhonov regularization using CGLS
        */        
        void regularize_tikhonov_CG(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(LQR, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
            dBlockType Tk(LQR, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

            // ----------------------------------------------------------------------------------------------
            // generate x0
            // ----------------------------------------------------------------------------------------------
            dVectorType sol_x(ObjRank+RemainingColumns);
            sol_x.setZero();
            // ----------------------------------------------------------------------------------------------

            cg_tikhonov(sol_x, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol_x.tail(RemainingColumns);
        }

        /** 
            \brief Tikhonov regularization using CGLS

            \note hot-start from RT_NO_Z
        */        
        void regularize_tikhonov_CG_x0(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            RealScalar mu = damp_factor*damp_factor;

            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(  LQR, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
            dBlockType Tk(  LQR, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);
            dBlockType  D(array,             0,                     0, ObjRank,          ObjRank);
            
            dBlockType2Vector d(array, 0, nVar, ObjRank, 1);

            // ----------------------------------------------------------------------------------------------
            // generate x0
            // ----------------------------------------------------------------------------------------------    
            D.triangularView<Eigen::Lower>() = (Rk.triangularView<Eigen::Upper>()*Rk.transpose()).eval();
            D.triangularView<Eigen::Lower>() += Tk*Tk.transpose();

            for (Index i=0; i<ObjRank; i++)
                D.coeffRef(i,i) += mu;

            d = LQR.col(nVar).segment(FirstRowIndex,ObjRank);

            Eigen::LLT<MatrixType> chol(D);
            chol.solveInPlace(d);

            dVectorType sol(ObjRank+RemainingColumns);
            sol.head(ObjRank).noalias() = Rk.triangularView<Eigen::Upper>().transpose() * d;
            sol.tail(RemainingColumns) = Tk.transpose() * d;
            // ----------------------------------------------------------------------------------------------

            cg_tikhonov(sol, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol.head(ObjRank);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol.tail(RemainingColumns);
        }

        /** 
            \brief RT_NO_Z regularization using CGLS
        */        
        void regularize_RT_NO_Z_CG(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            // -------------------------------------------------------------------------
            // create blocks
            // -------------------------------------------------------------------------
            dBlockType Rk(LQR, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
            dBlockType Tk(LQR, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

            // ----------------------------------------------------------------------------------------------
            // generate x0
            // ----------------------------------------------------------------------------------------------
            dVectorType sol(ObjRank+RemainingColumns);
            sol.setZero();
            // ----------------------------------------------------------------------------------------------

            cg_RT(sol, FirstRowIndex, FirstColIndex, ObjRank, RemainingColumns);

            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias()  = Rk.triangularView<Eigen::Upper>()*sol.head(ObjRank);
            LQR.col(nVar).segment(FirstRowIndex,ObjRank).noalias() += Tk*sol.tail(RemainingColumns);
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
            dBlockType Rk(      LQR, FirstRowIndex,         FirstColIndex,       ObjRank,                    ObjRank);
            dBlockType Tk(      LQR, FirstRowIndex, FirstColIndex+ObjRank,       ObjRank,           RemainingColumns);
            dBlockType Sk(NullSpace,             0,         FirstColIndex, FirstColIndex, RemainingColumns + ObjRank);

            dBlockType2Vector r(rqsp_work, 0, 0, ObjRank + FirstColIndex + ObjRank + RemainingColumns, 1);
            dBlockType2Vector q(rqsp_work, 0, 1, ObjRank + FirstColIndex + ObjRank + RemainingColumns, 1);
            dBlockType2Vector s(rqsp_work, 0, 2,                           ObjRank + RemainingColumns, 1);
            dBlockType2Vector p(rqsp_work, 0, 3,                           ObjRank + RemainingColumns, 1);
            
            // ------------------------------------------------------------------------------------------------
            /*
                  | yk |   | Rk Tk | 
              r = | sk | - |  Sk   | * x
                  |  0 |   |  Ik   |
            */
            // ------------------------------------------------------------------------------------------------
            r.head(ObjRank).noalias()  = LQR.col(nVar).segment(FirstRowIndex,ObjRank) - Tk*sol_x.tail(RemainingColumns);
            r.head(ObjRank).noalias() -= Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);

            r.segment(ObjRank,FirstColIndex).noalias() = NullSpace.col(nVar).head(FirstColIndex) - Sk * sol_x;
            r.segment(ObjRank,FirstColIndex)          *= damp_factor;

            r.tail(ObjRank+RemainingColumns).noalias() = -damp_factor*sol_x;
            // ------------------------------------------------------------------------------------------------
            /*
                  | Rk'       |   | r1 |
              s = |     Sk' I | * | r2 |
                  | Tk'       |   | r3 |
            */
            // ------------------------------------------------------------------------------------------------
            s.noalias() = Sk.transpose() * r.segment(ObjRank,FirstColIndex) + r.tail(ObjRank+RemainingColumns);
            s *= damp_factor;

            s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
            s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
            // ------------------------------------------------------------------------------------------------
            p = s;
            // ------------------------------------------------------------------------------------------------
            gamma = s.squaredNorm();
            // ------------------------------------------------------------------------------------------------

            Index iter = 0;
            while ( (std::sqrt(gamma) > tol) && (iter < regularizationMaxIterCG) )
            {
                // ------------------------------------------------------------------------------------------------
                // q = [Rk Tk; Sk ; Ik]*p
                // ------------------------------------------------------------------------------------------------
                q.head(ObjRank).noalias()  = Tk*p.tail(RemainingColumns);
                q.head(ObjRank).noalias() += Rk.triangularView<Eigen::Upper>()*p.head(ObjRank);

                q.segment(ObjRank,FirstColIndex).noalias() = Sk * p;
                q.segment(ObjRank,FirstColIndex) *= damp_factor;

                q.tail(ObjRank+RemainingColumns).noalias() = damp_factor*p;
                // ------------------------------------------------------------------------------------------------
                alpha = gamma/q.squaredNorm();
                sol_x += alpha*p;
                r -= alpha*q;
                // ------------------------------------------------------------------------------------------------
                // S = [Rk Tk; Sk ; Ik]'*r
                // ------------------------------------------------------------------------------------------------
                s.noalias() = Sk.transpose() * r.segment(ObjRank,FirstColIndex) + r.tail(ObjRank+RemainingColumns);
                s *= damp_factor;

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

            dBlockType Rk( LQR, FirstRowIndex,         FirstColIndex, ObjRank,          ObjRank);
            dBlockType Tk( LQR, FirstRowIndex, FirstColIndex+ObjRank, ObjRank, RemainingColumns);

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
            r.head(ObjRank).noalias()  = LQR.col(nVar).segment(FirstRowIndex,ObjRank) - Tk*sol_x.tail(RemainingColumns);
            r.head(ObjRank).noalias() -= Rk.triangularView<Eigen::Upper>()*sol_x.head(ObjRank);

            r.tail(ObjRank+RemainingColumns).noalias() = -damp_factor*sol_x;
            // ------------------------------------------------------------------------------------------------
            /*
                  | Rk'   |   | r1 |
              s = |     I | * | r2 |
                  | Tk'   |
            */
            // ------------------------------------------------------------------------------------------------
            s.noalias() = damp_factor * r.tail(ObjRank+RemainingColumns);

            s.head(ObjRank).noalias()          += Rk.triangularView<Eigen::Upper>().transpose() * r.head(ObjRank);
            s.tail(RemainingColumns).noalias() += Tk.transpose() * r.head(ObjRank);
            // ------------------------------------------------------------------------------------------------
            p = s;
            // ------------------------------------------------------------------------------------------------
            gamma = s.squaredNorm();
            // ------------------------------------------------------------------------------------------------

            Index iter = 0;
            while ( (std::sqrt(gamma) > tol) && (iter < regularizationMaxIterCG) )
            {
                // ------------------------------------------------------------------------------------------------
                // q = [Rk Tk; Ik]*p
                // ------------------------------------------------------------------------------------------------
                q.head(ObjRank).noalias()  = Tk*p.tail(RemainingColumns);
                q.head(ObjRank).noalias() += Rk.triangularView<Eigen::Upper>()*p.head(ObjRank);

                q.tail(ObjRank+RemainingColumns).noalias() = damp_factor*p;
                // ------------------------------------------------------------------------------------------------
                alpha = gamma/q.squaredNorm();
                sol_x += alpha*p;
                r -= alpha*q;
                // ------------------------------------------------------------------------------------------------
                // S = [Rk Tk; Ik]'*r
                // ------------------------------------------------------------------------------------------------
                s.noalias() = damp_factor * r.tail(ObjRank+RemainingColumns);

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
            \brief Accumulate nullspace basis
        */
        void accumulate_nullspace_basis(Index FirstRowIndex, Index FirstColIndex, Index ObjRank, Index RemainingColumns)
        {
            dBlockType            Rk(      LQR, FirstRowIndex,           FirstColIndex,               ObjRank,            ObjRank);
            dBlockType       UpBlock(      LQR, FirstRowIndex, FirstColIndex + ObjRank,               ObjRank, RemainingColumns+1);
            dBlockType     LeftBlock(NullSpace,             0,           FirstColIndex, FirstColIndex+ObjRank,            ObjRank);
            dBlockType TrailingBlock(NullSpace,             0, FirstColIndex + ObjRank, FirstColIndex+ObjRank, RemainingColumns+1);

            // I have to do this because in regularize_tikhonov_1(...) I use NullSpace as a temporary storage
            LeftBlock.block(FirstColIndex,0,ObjRank,ObjRank).setIdentity();

            Rk.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(LeftBlock);
            
            if (ObjRank == 1)
            {
                TrailingBlock.noalias() -= LeftBlock.col(0) * UpBlock.row(0);
            }
            else if (ObjRank >= 2 && ObjRank <= 8)
            {
                for (Index k=0; k<RemainingColumns+1; k++)
                    TrailingBlock.col(k).noalias() -= LeftBlock * UpBlock.col(k);
            }                        
            else if (ObjRank > 8)
            {
                TrailingBlock.noalias() -= LeftBlock * UpBlock;
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
            \brief Max number of iterations for cg_tikhonov(...)

            \note used only with regularizationType = REGULARIZATION_TIKHONOV_CG
        */
        Index regularizationMaxIterCG;
  
        /** 
            \brief Linear dependence tolerance
        */
        RealScalar LinearDependenceTolerance;

        /** 
            \brief Regularization factor used at the current level

            \note determined as a function of regularization[ObjIndex] and the smallest pivot
        */
        RealScalar damp_factor;

        // ==================================================================
        // bool
        // ==================================================================

        /** 
            \brief isFactorized = false while the factorization has not been completed
        */
        bool isFactorized;

        /** 
            \brief isSolved = false while the problem has not been solved
        */
        bool isSolved;

        /** 
            \brief use the real residual when computing sensitivity (or not)
        */
        bool realSensitivityResidual;

        // ==================================================================
        // definition of vectors of integers
        // ==================================================================

        /** 
            \brief Stores the indexes of permuted columns
        */
        iVectorType ColPermutations;

        /** 
            \brief Indexes of fixed variables

            \attention This field is changed in #factorize().
        */
        iVectorType FixedVarIndex;

        /** 
            \brief For debugging purposes only (FIXME: to remove) 

            \note Cannot be used with simpli-bounded first objective
        */
        iVectorType DedicatedVariables;

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

            \note #LQR and #x are kept separately for convenience. 
        */
        dVectorType dWorkspace;

        /** 
            \brief A (primal) solution

            \note When TotalRank < nVar the problem is under-determined and x is a "basic solution" 
            (see p. 106 of "Numerical Methods for Least Squares Problems" by Bjorck, 1996).
        */
        dVectorType x;

        /*
          \note Householder scalars associated with the QR factorization for a (LexLSE) objective. 
          
          \note computed during the factorization.
        */
        dVectorType hh_scalars; 

        /** 
            \brief A heuristic regularization for each objective

            \note The number of elements is equal to the number of objectives (fixed variables in
            the highest level are never regularized)
        */
        dVectorType regularization; 

        /** 
            \brief Used in implementation of conjugate gradients (CG)
        */
        MatrixType rqsp_work;
        
        // ==================================================================
        // definition of matrices
        // ==================================================================

        /** 
            \brief Constraint matrix 

            \note On input, LQR contains the matrix of constraints (the last column contains the RHS
            vector) of the lexicographic problem. On output, LQR contains the matrices L, Q and R
            from the LQR factorization of the constraint matrix. The orthonormal matrices Q are
            stored in Householder factored form below the main diagonal of the leading triangular
            matrices of each objective, while the last column contains (LQ)^{-1}*RHS.

            \note The row-size of LQR is set to the the maximum number of envisioned constraints
            (excluding "fixing constraints"). This is usefull in the context of LexLSI where one
            instance of LexLSE is used to solve problems with different dimensoins.
        */
        MatrixType LQR;

        /** 
            \brief nVar x (nVar+1) matrix containing the remaining null-space + right-hand-side

            \note Used for Tikhonov regularization or terminal objective of the form (x = 0)
        */        
        MatrixType NullSpace;

        /** 
            \brief Used for temporary storage when doing Tikhonov regularization (option: A'*inv(A*A'+I)*b)
        */        
        MatrixType array;

        // ==================================================================
        // other
        // ==================================================================

        /** 
            \brief Vector containing information about the constraints involved in each objective
        */
        std::vector<ObjInfoType> ObjInfo;

        /** 
            \brief Type of fixed variables (LOWER_BOUND or UPPER_BOUND)
        */
        std::vector<ConstraintType> FixedVarType;

        /** 
            \brief Type of (active) constraints (LOWER_BOUND or UPPER_BOUND)
        */
        std::vector<ConstraintType> CtrType;

        /** 
            \brief Permutation matrix
        */
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;

        /** 
            \brief Type of regularization (Tikhonov, Basic Tikhonov, ...)
        */
        RegularizationType regularizationType;

    }; // end class LexLSE

} // end namespace Eigen

#endif // LEXLSE
