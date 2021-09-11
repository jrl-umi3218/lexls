/*
 * Copyright 2013-2021 INRIA
 */

/**
    @file
    @author  Alexander Sherikov

    @brief 
*/


#include <lexlsi.h>
#include <tools.h>

#include <iostream>
#include <iomanip>


int main(int argc, char **argv)
{
    if (argc != 2)
    {
        return (-1);
    }
    

    LexLS::tools::HierarchyType          type_of_hierarchy;
    unsigned int                         number_of_variables;
    unsigned int                         number_of_objectives;
    std::vector<unsigned int>            number_of_constraints;
    std::vector<LexLS::ObjectiveType>    types_of_objectives;
    std::vector<Eigen::MatrixXd>         objectives;

    std::vector< std::vector<LexLS::ConstraintActivationType> >  active_set_guess;
    Eigen::VectorXd                                     solution_guess;
    Eigen::VectorXd                                     solution;

    try
    {
        LexLS::tools::HierarchyFileProcessor fprocessor;

        fprocessor.import(  argv[1],
                            type_of_hierarchy,
                            number_of_variables,
                            number_of_objectives,
                            number_of_constraints,
                            types_of_objectives,
                            objectives,
                            active_set_guess,
                            solution_guess,
                            solution);

        // ==============================================================
        for (unsigned int i = 0; i < number_of_objectives; ++i)
        {
            std::cout << objectives[i] << std::endl << std::endl << std::endl << std::endl;
        }
        // ==============================================================


        // ==============================================================
        if (active_set_guess.size() != 0)
        {
            for (unsigned int i = 0; i < number_of_objectives; ++i)
            {
                for (unsigned int j = 0; j < number_of_constraints[i]; ++j)
                {
                    std::cout << active_set_guess[i][j] << std::endl;
                }
                std::cout << std::endl << std::endl << std::endl;
            }
        }
        // ==============================================================

    
        // ==============================================================
        if (solution_guess.size() != 0)
        {
            for (unsigned int i = 0; i < number_of_variables; ++i)
            {
                std::cout << solution_guess[i] << std::endl;
            }
            std::cout << std::endl << std::endl << std::endl;
        }
        // ==============================================================

        
        // ==============================================================
        if (solution.size() != 0)
        {
            for (unsigned int i = 0; i < number_of_variables; ++i)
            {
                std::cout << solution[i] << std::endl;
            }
            std::cout << std::endl << std::endl << std::endl;
        }
        // ==============================================================
    }
    catch(const std::exception &e)
    {
        std::cout << e.what() << std::endl;
    }


    return (0);
}
