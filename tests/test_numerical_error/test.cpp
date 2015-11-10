/*
  Problem:
  ---------
  In general the solution of the last equality constrained problem (solved in lexlsi) cannot be
  reproduced by passing exactly the same data to lexlse manually. It seems that the compiler treats
  the code differently (maybe some computations are sequenced differently). As a result there are
  small errors. When we have many levels (as in Nestor's examples) this error magnifies (I have an
  example where the norm of the difference is greater than 1e-02).

  Solution:
  ---------
  This is due to Eigen's vectorization. If one defined EIGEN_DONT_VECTORIZE, the problem is gone
  (but we don't want to do that).
*/

#include <lexlsi.h>
#include <tools.h>

#include <iostream>
#include <iomanip>

int main()
{

    LexLS::tools::HierarchyType       type_of_hierarchy;
    LexLS::Index                      number_of_variables;
    LexLS::Index                      number_of_objectives;
    std::vector<LexLS::Index>         number_of_constraints;
    std::vector<LexLS::ObjectiveType> types_of_objectives;
    std::vector<Eigen::MatrixXd>      objectives;

    std::vector< std::vector<LexLS::ConstraintActivationType> > active_set_guess;
    Eigen::VectorXd solution_guess;
    Eigen::VectorXd solution;

    // ==============================================================

    LexLS::tools::HierarchyFileProcessor fprocessor;
    fprocessor.import("test_113.dat",
                      type_of_hierarchy,
                      number_of_variables,
                      number_of_objectives,
                      number_of_constraints,
                      types_of_objectives,
                      objectives,
                      active_set_guess,
                      solution_guess,
                      solution);

    LexLS::internal::LexLSI lsi(number_of_variables,
                                number_of_objectives,
                                &number_of_constraints[0],
                                &types_of_objectives[0]);

    LexLS::ParametersLexLSI parameters;
    parameters.max_number_of_factorizations = 113;
    lsi.setParameters(parameters);

    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        lsi.setData(i,objectives[i]);
    }

    lsi.solve();

    int nCtr = 0;
    std::vector<LexLS::Index> eq_obj_dim(number_of_objectives);
    for (LexLS::Index k=0; k<number_of_objectives; k++)
    {
        eq_obj_dim[k] = lsi.getActiveCtrCount(k);

        nCtr += eq_obj_dim[k];
        std::cout << eq_obj_dim[k] << " ";
    }
    std::cout << "\n\n";

    // ==============================================================

    LexLS::dMatrixType X(number_of_variables,3);

    X.col(0) = lsi.get_xStar();

    LexLS::dMatrixType data = lsi.get_data(); // call after lsi.get_xStar()
    LexLS::dMatrixType lexqr = lsi.get_lexqr();

    std::vector<LexLS::ConstraintIdentifier> ctr;
    lsi.getActiveCtr_order(ctr);

    // ==============================================================
    LexLS::internal::LexLSE lse(number_of_variables,
                                number_of_objectives,
                                &eq_obj_dim[0]);

    LexLS::Index row_ind = 0;
    int k=0;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        lse.setData(i,data.block(row_ind, 0, eq_obj_dim[i], number_of_variables+1));

        for (unsigned int j = 0; j < eq_obj_dim[i]; ++j)
        {
            lse.setCtrType(i, j, ctr[k].ctr_type);
            k++;
        }

        row_ind += eq_obj_dim[i];
    }

    // ==============================================================
    lse.factorize();
    lse.solve();

    // ==============================================================
    LexLS::dMatrixType data1 = lse.get_data();
    LexLS::dMatrixType lexqr1 = lse.get_lexqr();

    //std::cout << "nCtr = " << nCtr << std::endl;
    //std::cout << data.rows() << ", " << data.cols() << std::endl;
    //std::cout << data1.rows() << ", " << data1.cols() << std::endl;

    LexLS::RealScalar e1 = (data.topRows(nCtr) - data1).norm();

    LexLS::RealScalar e2 = (lexqr.topRows(nCtr) - lexqr1).norm();
    // ==============================================================

    X.col(1) = lse.get_x();
    X.col(2) = X.col(1) - X.col(0);

    //std::cout << X << "\n\n";

    printf("error(x)     = %1.18e \n", X.col(2).norm());
    printf("error(data)  = %1.18e \n", e1);
    printf("error(lexqr) = %1.18e \n", e2);

    // ==============================================================
    // residual
    // ==============================================================
/*
    LexLS::dMatrixType R(nCtr,3);

    std::cout << data.rows() << ", " << data.cols() << std::endl;
    std::cout << number_of_variables << std::endl;
    std::cout << X.col(1).size() << std::endl;

    R.col(0) = data.topRows(nCtr).leftCols(number_of_variables)*X.col(0) - data.topRows(nCtr).col(number_of_variables);
    R.col(1) = data1.leftCols(number_of_variables)*X.col(1) - data1.col(number_of_variables);
    R.col(2) = R.col(1) - R.col(0);


    row_ind = 0;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        std::cout << "----------------- obj = " << i << "-----------------\n";
        std::cout << R.block(row_ind, 0, eq_obj_dim[i], 3) << "\n";

        row_ind += eq_obj_dim[i];
    }
*/

    row_ind = 0;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        double e = (lexqr.block(row_ind, 0, eq_obj_dim[i], number_of_variables+1) - lexqr1.block(row_ind, 0, eq_obj_dim[i], number_of_variables+1)).norm();
        printf("   obj(%d): %1.18e \n",i,e);

        row_ind += eq_obj_dim[i];
    }

    return 0;
}
