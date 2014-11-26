// Time-stamp: <2014-11-26 11:07:20 drdv>
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
            
            Index dim = std::max(maxObjDimSum,nVar);
            dWorkspace.resize(2*dim + nVar + 1);
            
            regularization.resize(nObj);
            ColPermutations.resize(nVar);
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

                // -----------------------------------------------------------------------
                // Gauss transformation
                // -----------------------------------------------------------------------
                ObjRank = ObjInfo[ObjIndex].rank = ColIndex - FirstColIndex; // store the rank

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

                        // ====================================================================================================
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
                                  RealScalar tolWrongSignLambda, RealScalar tolCorrectSignLambda)
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

            P.setIdentity();
            x.tail(nVar-nVarFixed).setZero(); // x.head(nVarFixed) has already been initialized in fixVariable(...)

            regularization.setZero(); // by default there is no regularization
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
            \brief Linear dependence tolerance
        */
        RealScalar LinearDependenceTolerance;

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

    }; // end class LexLSE

} // end namespace Eigen

#endif // LEXLSE
