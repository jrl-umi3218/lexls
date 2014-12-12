/**
    @file
    @author  Alexander Sherikov

    @brief 
*/


#include <lexlsi.h>
#include <lexls_tools.h>

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
    std::vector<unsigned int>            number_of_constraits_per_objective;
    std::vector<LexLS::ObjectiveType>    types_of_objectives;
    std::vector<Eigen::MatrixXd>         objectives;


    try
    {
        import_hierarchy(  argv[1],
                           type_of_hierarchy,
                           number_of_variables,
                           number_of_objectives,
                           number_of_constraits_per_objective,
                           types_of_objectives,
                           objectives);

        // ==============================================================
        for (unsigned int i = 0; i < number_of_objectives; ++i)
        {
            std::cout << objectives[i] << std::endl << std::endl << std::endl << std::endl;
        }
        // ==============================================================
    }
    catch(const std::exception &e)
    {
        std::cout << e.what() << std::endl;
    }


    return (0);
}
