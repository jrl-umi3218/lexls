/*
  In general the solution of the last equality constrained problem (solved in lexlsi) cannot be
  reproduced by passing exactly the same data to lexlse manually. It seems that the compiler treats
  the code differently (maybe some computations are sequenced differently). As a result there are
  small errors. When we have many levels (as in Nestor's examples) this error magnifies (I have an
  example where the norm of the difference is greater than 1e-02).
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

    std::vector<LexLS::Index> eq_obj_dim(number_of_objectives);
    for (LexLS::Index k=0; k<number_of_objectives; k++)
    {
        eq_obj_dim[k] = lsi.getActiveCtrCount(k);

        std::cout << eq_obj_dim[k] << " ";
    }
    std::cout << "\n\n";

    // ==============================================================

    LexLS::dMatrixType X(number_of_variables,3);

    X.col(0) = lsi.get_xStar();

    LexLS::dMatrixType data = lsi.get_data(); // call after lsi.get_xStar()

    // ==============================================================
    LexLS::internal::LexLSE lse(number_of_variables,
                                number_of_objectives,
                                &eq_obj_dim[0]);

    LexLS::Index row_ind = 0;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        if (i > 0)
        {
            row_ind += eq_obj_dim[i-1];
        }

        lse.setData(i,data.block(row_ind, 0, eq_obj_dim[i], number_of_variables+1));
    }
    // ==============================================================
    lse.factorize();
    lse.solve();

    X.col(1) = lse.get_x();
    X.col(2) = X.col(1) - X.col(0);

    std::cout << X << "\n\n";

    printf("% 1.18e \n", X.col(2).norm());

    return 0;
}
