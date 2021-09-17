/*
 * Copyright 2013-2021 INRIA
 */

#pragma once

#include <Eigen/Core>
#include <stdexcept>
#include <vector>
#include <fstream>

namespace LexLS
{
    namespace tools
    {
        enum HierarchyType
        {
            HIERARCHY_TYPE_NONE         = 0,
            HIERARCHY_TYPE_EQUALITY     = 1, // equality constraints
            HIERARCHY_TYPE_INEQUALITY   = 2 // inequality constraints
        };


        class HierarchyFileProcessor
        {
            private:
                enum HierarchyTypeHeader
                {
                    HEADER_HIERARCHY_NONE                  = 0,   //
                    HEADER_HIERARCHY_EQUALITIES            = 100, // equality constraints
                    HEADER_HIERARCHY_INEQUALITIES          = 200, // inequality constraints
                    HEADER_HIERARCHY_INEQUALITIES_WITH_AS  = 210  // inequality constraints with active set guess
                };


                enum ObjectiveTypeHeader
                {
                    HEADER_OBJECTIVE_TYPE_SIMPLE   = 100,
                    HEADER_OBJECTIVE_TYPE_GENERAL  = 200
                };


                std::string NUMBER_OF_VARIABLES  ;
                std::string NUMBER_OF_OBJECTIVES ;
                std::string NUMBER_OF_CONSTRAINTS;
                std::string HIERARCHY_TYPE       ;
                std::string TYPES_OF_OBJECTIVES  ;
                std::string OBJECTIVE_DATA       ;
                std::string SOLUTION_GUESS       ;
                std::string SOLUTION             ;


                void readNumberOfVariables(
                        std::ifstream   & ifs,
                        bool            & set_number_of_variables,
                        unsigned int    & number_of_variables)
                {
                    if (set_number_of_variables)
                    {
                        throw std::runtime_error("Duplicate header field.");
                    }

                    ifs >> number_of_variables;

                    if (!ifs.fail())
                    {
                        set_number_of_variables = true;
                    }
                }


                void readNumberOfObjectives(
                        std::ifstream   & ifs,
                        bool            & set_number_of_objectives,
                        unsigned int    & number_of_objectives)
                {
                    if (set_number_of_objectives)
                    {
                        throw std::runtime_error("Duplicate header field.");
                    }

                    ifs >> number_of_objectives;
                    if (!ifs.fail())
                    {
                        set_number_of_objectives = true;
                    }
                }


                void readTypeOfHieararchy(
                        std::ifstream   & ifs,
                        bool            & set_type_of_hierarchy,
                        HierarchyType   & type_of_hierarchy,
                        unsigned int    & type_of_hierarchy_header)
                {
                    if (set_type_of_hierarchy)
                    {
                        throw std::runtime_error("Duplicate header field.");
                    }


                    ifs >> type_of_hierarchy_header;


                    if (!ifs.fail())
                    {
                        set_type_of_hierarchy = true;
                        switch(type_of_hierarchy_header)
                        {
                            case HEADER_HIERARCHY_EQUALITIES:
                                type_of_hierarchy = HIERARCHY_TYPE_EQUALITY;
                                break;
                            case HEADER_HIERARCHY_INEQUALITIES:
                                type_of_hierarchy = HIERARCHY_TYPE_INEQUALITY;
                                break;
                            case HEADER_HIERARCHY_INEQUALITIES_WITH_AS:
                                type_of_hierarchy = HIERARCHY_TYPE_INEQUALITY;
                                break;
                            default:
                                throw std::runtime_error("Unsupported type of hierarchy.");
                                break;
                        }
                    }
                }


                void readNumberOfConstraints(
                        std::ifstream               & ifs,
                        bool                        & set_number_of_constraints,
                        std::vector<unsigned int>   & number_of_constraints)
                {
                    if (set_number_of_constraints)
                    {
                        throw std::runtime_error("Duplicate header field.");
                    }


                    std::string     line; // current line
                    getline(ifs, line);
                    std::stringstream   line_stream;
                    line_stream.str(line);
                    line_stream.clear();

                    while (!line_stream.eof())
                    {
                        unsigned int ctr_num = 0;

                        line_stream >> ctr_num;

                        if (!line_stream.fail())
                        {
                            number_of_constraints.push_back(ctr_num);
                        }
                    }

                    set_number_of_constraints = true;
                }


                void readTypesOfObjectives(
                        std::ifstream                       & ifs,
                        bool                                & set_types_of_objectives,
                        std::vector<LexLS::ObjectiveType>   & types_of_objectives)
                {
                    if (set_types_of_objectives)
                    {
                        throw std::runtime_error("Duplicate header field.");
                    }


                    std::string     line; // current line
                    getline(ifs, line);
                    std::stringstream   line_stream;
                    line_stream.str(line);
                    line_stream.clear();

                    while (!line_stream.eof())
                    {
                        unsigned int type;

                        line_stream >> type;

                        if (!line_stream.fail())
                        {
                            switch(type)
                            {
                                case HEADER_OBJECTIVE_TYPE_SIMPLE:
                                    types_of_objectives.push_back(LexLS::SIMPLE_BOUNDS_OBJECTIVE);
                                    break;
                                case HEADER_OBJECTIVE_TYPE_GENERAL:
                                    types_of_objectives.push_back(LexLS::GENERAL_OBJECTIVE);
                                    break;
                                default:
                                    throw std::runtime_error("Unsupported type of objective.");
                                    break;
                            }
                        }
                    }

                    set_types_of_objectives = true;
                }


                void readObjective( const unsigned int  number_of_rows,
                                    const unsigned int  number_of_columns,
                                    const unsigned int  objective_index,
                                    const unsigned int  type_of_hierarchy_header,
                                    std::ifstream               & ifs,
                                    std::vector<Eigen::MatrixXd>                        & objectives,
                                    std::vector< std::vector<ConstraintActivationType> >  & active_set_guess)
                {
                    // excessive rows / columns are silently ignored.
                    for (unsigned int row_index = 0;  row_index < number_of_rows;  ++row_index)
                    {
                        std::string     line; // current line

                        getline(ifs, line);
                        std::stringstream   line_stream;
                        line_stream.str(line);
                        line_stream.clear();

                        unsigned int column_index = 0;
                        while ( !line_stream.eof()  &&
                                !line_stream.fail()      &&
                                (column_index < number_of_columns))
                        {
                            line_stream >> objectives[objective_index](row_index, column_index);
                            ++column_index;
                        }


                        if (column_index != number_of_columns)
                        {
                            throw std::runtime_error("Not enough data.");
                        }


                        if (type_of_hierarchy_header == HEADER_HIERARCHY_INEQUALITIES_WITH_AS)
                        {
                            unsigned int    activation_type;
                            line_stream >>  activation_type;

                            if (!line_stream.fail())
                            {
                                switch(activation_type)
                                {
                                    case CTR_INACTIVE:
                                    case CTR_ACTIVE_LB:
                                    case CTR_ACTIVE_UB:
                                    case CTR_ACTIVE_EQ:
                                        active_set_guess[objective_index][row_index] =
                                            static_cast<ConstraintActivationType> (activation_type);
                                        break;
                                    default:
                                        throw std::runtime_error("Unsupported constraint activation type.");
                                        break;
                                }
                            }
                        }
                    }
                }


                void readSolution(  const unsigned int  number_of_variables,
                                    std::ifstream       & ifs,
                                    Eigen::VectorXd     & solution_vector)
                {
                    solution_vector.resize(number_of_variables);

                    for (unsigned int i = 0; i < number_of_variables; ++i)
                    {
                        ifs >> solution_vector(i);

                        if (ifs.fail())
                        {
                            throw std::runtime_error("Could not read a solution vector.");
                        }
                    }
                }



            public:
                HierarchyFileProcessor()
                {
                    NUMBER_OF_VARIABLES         = "#nVar";
                    NUMBER_OF_OBJECTIVES        = "#nObj";
                    NUMBER_OF_CONSTRAINTS       = "#nCtr";
                    HIERARCHY_TYPE              = "#HierType";
                    TYPES_OF_OBJECTIVES         = "#ObjType";
                    OBJECTIVE_DATA              = "#OBJECTIVE";
                    SOLUTION_GUESS              = "#SolGuess";
                    SOLUTION                    = "#Solution";
                }


                void import(const std::string                       & file_name,
                            HierarchyType                           & type_of_hierarchy,
                            unsigned int                            & number_of_variables,
                            unsigned int                            & number_of_objectives,
                            std::vector<unsigned int>               & number_of_constraints,
                            std::vector<LexLS::ObjectiveType>       & types_of_objectives,
                            std::vector<Eigen::MatrixXd>            & objectives,
                            std::vector< std::vector<ConstraintActivationType> >   & active_set_guess,
                            Eigen::VectorXd                         & solution_guess,
                            Eigen::VectorXd                         & solution)
                {
                    type_of_hierarchy = HIERARCHY_TYPE_NONE;
                    number_of_variables = 0;
                    number_of_objectives = 0;

                    number_of_constraints.resize(0);
                    types_of_objectives.resize(0);
                    objectives.resize(0);
                    active_set_guess.resize(0);
                    solution_guess.resize(0);
                    solution.resize(0);


                    std::ifstream ifs;
                    ifs.open(file_name.c_str());

                    if (!ifs)
                    {
                        throw std::runtime_error("Cannot open file for reading");
                    }


                    // ---------------------------------------------------------------------------


                    bool set_number_of_variables    = false;
                    bool set_number_of_objectives   = false;
                    bool set_number_of_constraints  = false;
                    bool set_type_of_hierarchy      = false;
                    bool set_types_of_objectives    = false;

                    unsigned int type_of_hierarchy_header = HEADER_HIERARCHY_NONE;

                    // ---------------------------------------------------------------------------

                    std::string     line; // current line

                    // read header, fields in the header can be in an arbitrary order
                    while (getline(ifs, line) &&
                            ! ( set_number_of_variables && set_number_of_objectives &&
                                set_type_of_hierarchy   && set_types_of_objectives &&
                                set_number_of_constraints) )
                    {
                        // remove empty space from a string (just in case)
                        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

                        if (line.compare(NUMBER_OF_VARIABLES) == 0)
                        {
                            readNumberOfVariables(ifs, set_number_of_variables, number_of_variables);
                            continue;
                        }

                        if (line.compare(NUMBER_OF_OBJECTIVES) == 0)
                        {
                            readNumberOfObjectives(ifs, set_number_of_objectives, number_of_objectives);
                            continue;
                        }

                        if (line.compare(HIERARCHY_TYPE) == 0)
                        {
                            readTypeOfHieararchy(ifs, set_type_of_hierarchy, type_of_hierarchy, type_of_hierarchy_header);
                            continue;
                        }

                        if (line.compare(NUMBER_OF_CONSTRAINTS) == 0)
                        {
                            readNumberOfConstraints(ifs, set_number_of_constraints, number_of_constraints);
                            continue;
                        }

                        if (line.compare(TYPES_OF_OBJECTIVES) == 0)
                        {
                            readTypesOfObjectives(ifs, set_types_of_objectives, types_of_objectives);
                            continue;
                        }
                    }


                    if (    !set_number_of_variables        || !set_number_of_objectives ||
                            !set_type_of_hierarchy          || !set_types_of_objectives ||
                            !set_number_of_constraints)
                    {
                        throw std::runtime_error("At least one required parameters is not set.");
                    }


                    if (    (types_of_objectives.size() != number_of_objectives) ||
                            (number_of_constraints.size() != number_of_objectives))
                    {
                        throw std::runtime_error("Wrong number of objectives.");
                    }


                    // ---------------------------------------------------------------------------


                    objectives.resize(number_of_objectives);

                    // find OBJECTIVE_DATA
                    // (objectives are assumed to be stored in ascending order after header fields)
                    unsigned int  objective_index = 0;
                    unsigned int  number_of_columns = 0;
                    unsigned int  number_of_bounds = 0;

                    switch (type_of_hierarchy_header)
                    {
                        case HEADER_HIERARCHY_EQUALITIES:
                            number_of_bounds = 1;
                            break;
                        case HEADER_HIERARCHY_INEQUALITIES:
                            number_of_bounds = 2;
                            break;
                        case HEADER_HIERARCHY_INEQUALITIES_WITH_AS:
                            active_set_guess.resize(number_of_objectives);
                            number_of_bounds = 2;
                            break;
                        default:
                            throw std::runtime_error("Unsupported type of hierarchy.");
                            break;
                    }


                    while (getline(ifs,line) && (objective_index < number_of_objectives))
                    {
                        // remove empty space from a string (just in case)
                        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

                        if (line.compare(0, OBJECTIVE_DATA.length(), OBJECTIVE_DATA) == 0)
                        {
                            if ( types_of_objectives[objective_index] == LexLS::SIMPLE_BOUNDS_OBJECTIVE )
                            {
                                if (objective_index == 0)
                                {
                                    number_of_columns = 1 + number_of_bounds;
                                }
                                else
                                {
                                    throw std::runtime_error("Simple constraints are supported only in the first objective.");
                                }
                            }
                            else
                            {
                                number_of_columns = number_of_variables + number_of_bounds;
                            }

                            objectives[objective_index].resize(number_of_constraints[objective_index], number_of_columns);
                            if (type_of_hierarchy_header == HEADER_HIERARCHY_INEQUALITIES_WITH_AS)
                            {
                                active_set_guess[objective_index].resize(number_of_constraints[objective_index]);
                            }

                            readObjective(  number_of_constraints[objective_index],
                                            number_of_columns,
                                            objective_index,
                                            type_of_hierarchy_header,
                                            ifs,
                                            objectives,
                                            active_set_guess);

                            ++objective_index;
                        }
                    }


                    if (objective_index != number_of_objectives)
                    {
                        throw std::runtime_error("The number of objectives is lower than expected.");
                    }


                    while (getline(ifs, line))
                    {
                        if (line.compare(SOLUTION_GUESS) == 0)
                        {
                            readSolution(number_of_variables, ifs, solution_guess);
                            continue;
                        }

                        if (line.compare(SOLUTION) == 0)
                        {
                            readSolution(number_of_variables, ifs, solution);
                            continue;
                        }
                    }


                    ifs.close();
                }
        };

    } // END namespace internal

} // END namespace LexLS
