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
            HIERARCHY_EQUALITIES    = 1,
            HIERARCHY_INEQUALITIES  = 2
        };


        void import_hierarchy(  const std::string                   file_name,
                                HierarchyType                       & type_of_hierarchy,
                                unsigned int                        & number_of_variables,
                                unsigned int                        & number_of_objectives,
                                std::vector<unsigned int>           & number_of_constraits_per_objective,
                                std::vector<LexLS::ObjectiveType>   & types_of_objectives,
                                std::vector<Eigen::MatrixXd>        & objectives)
        {
            // Define identifiers
            const std::string NVAR = "#nVar";
            const std::string NOBJ = "#nObj";
            const std::string NCTR = "#nCtr";
            const std::string TYPE = "#Type";
            const std::string OBJTYPE = "#ObjType";
            const std::string DATA = "#OBJECTIVE";


            const unsigned int  OBJECTIVE_TYPE_SIMPLE   = 100;
            const unsigned int  OBJECTIVE_TYPE_GENERAL  = 200;
            

            // ---------------------------------------------------------------------------


            std::ifstream ifs;
            ifs.open(file_name.c_str());

            if (!ifs)
            {
                throw std::runtime_error("Cannot open file for reading");
            }


            // ---------------------------------------------------------------------------


            bool set_number_of_variables                = false;
            bool set_number_of_objectives               = false;
            bool set_number_of_constraits_per_objective = false;
            bool set_type_of_hierarchy                  = false;
            bool set_types_of_objectives                = false;


            // ---------------------------------------------------------------------------

            std::string     line; // current line

            // read header, fields in the header can be in an arbitrary order
            while (getline(ifs, line))
            {
                if (    set_number_of_variables &&
                        set_number_of_objectives && 
                        set_type_of_hierarchy && 
                        set_types_of_objectives && 
                        set_number_of_constraits_per_objective)
                {
                    break;
                }

                // remove empty space from a string (just in case)
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());


                if (line.compare(NVAR) == 0)
                {
                    if (!set_number_of_variables)
                    {
                        set_number_of_variables = true;
                        ifs >> number_of_variables;
                    }
                    else
                    {
                        throw std::runtime_error("Duplicated header field.");
                    }

                    continue;
                }
                

                if (line.compare(NOBJ) == 0)
                {
                    if (!set_number_of_objectives)
                    {
                        set_number_of_objectives = true;
                        ifs >> number_of_objectives;
                    }
                    else
                    {
                        throw std::runtime_error("Duplicated header field.");
                    }

                    continue;
                }


                if (line.compare(TYPE) == 0)
                {
                    if (!set_type_of_hierarchy)
                    {
                        set_type_of_hierarchy = true;

                        unsigned int type_of_hierarchy_tmp = 0;
                        ifs >> type_of_hierarchy_tmp;


                        if (    (type_of_hierarchy_tmp != HIERARCHY_EQUALITIES) && 
                                (type_of_hierarchy_tmp != HIERARCHY_INEQUALITIES))
                        {
                            throw std::runtime_error("Unsupported type of hierarchy.");
                        }

                        type_of_hierarchy = static_cast <HierarchyType> (type_of_hierarchy_tmp);
                    }
                    else
                    {
                        throw std::runtime_error("Duplicated header field.");
                    }

                    continue;
                }


                if (line.compare(NCTR) == 0)
                {
                    if (!set_number_of_constraits_per_objective)
                    {
                        std::stringstream   stream;
                        
                        getline(ifs, line);
                        stream.str(line);
                        stream.clear();


                        while (!stream.eof())
                        {
                            unsigned int ctr_num = 0;

                            stream >> ctr_num;

                            if (!stream.fail())
                            {
                                number_of_constraits_per_objective.push_back(ctr_num);
                            }
                        }
                        
                        set_number_of_constraits_per_objective = true;
                    }
                    else
                    {
                        throw std::runtime_error("Duplicated header field.");
                    }

                    continue;
                }


                if (line.compare(OBJTYPE) == 0)
                {
                    if (!set_types_of_objectives)
                    {
                        std::stringstream   stream;
                        
                        getline(ifs, line);
                        stream.str(line);
                        stream.clear();


                        while (!stream.eof())
                        {
                            unsigned int type;

                            stream >> type;

                            if (!stream.fail())
                            {
                                switch(type)
                                {
                                    case OBJECTIVE_TYPE_SIMPLE:
                                        types_of_objectives.push_back(LexLS::SIMPLE_BOUNDS_OBJECTIVE_HP);
                                        break;
                                    case OBJECTIVE_TYPE_GENERAL:
                                        types_of_objectives.push_back(LexLS::DEFAULT_OBJECTIVE);
                                        break;
                                    default:
                                        throw std::runtime_error("Unsupported type of objective.");
                                        break;
                                }
                            }
                        }

                        set_types_of_objectives = true;
                    }
                    else
                    {
                        throw std::runtime_error("Duplicated header field.");
                    }

                    continue;
                }
            }


            if (    !set_number_of_variables || 
                    !set_number_of_objectives || 
                    !set_type_of_hierarchy || 
                    !set_types_of_objectives || 
                    !set_number_of_constraits_per_objective)
            {
                throw std::runtime_error("At least one required parameters is not set.");
            }


            if (    (types_of_objectives.size() != number_of_objectives) || 
                    (number_of_constraits_per_objective.size() != number_of_objectives))
            {
                throw std::runtime_error("Wrong number of objectives.");
            }
           

            // ---------------------------------------------------------------------------


            objectives.resize(number_of_objectives);

            // find DATA (objectives are assumed to be stored in ascending order after header fields)
            unsigned int  objective_index = 0;
            unsigned int  number_of_columns = 0;
            unsigned int  number_of_bounds = 0;
            
            switch (type_of_hierarchy)
            {
                case HIERARCHY_EQUALITIES:
                    number_of_bounds = 1;
                    break;
                case HIERARCHY_INEQUALITIES:
                    number_of_bounds = 2;
                    break;
                default:
                    throw std::runtime_error("Unsupported type of hierarchy.");
                    break;
            }


            while (getline(ifs,line))
            {
                // remove empty space from a string (just in case)
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

                if (line.compare(0, DATA.length(), DATA) == 0)
                {
                    if ((objective_index != 0) || (types_of_objectives[objective_index] != LexLS::SIMPLE_BOUNDS_OBJECTIVE_HP))
                    {
                        number_of_columns = number_of_variables + number_of_bounds;
                    }
                    else
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

                    objectives[objective_index].resize(number_of_constraits_per_objective[objective_index], number_of_columns);

                    // excessive rows / columns are silently ignored.
                    for (unsigned int row_index = 0; row_index < number_of_constraits_per_objective[objective_index]; ++row_index)
                    {
                        std::stringstream   stream;
                        
                        getline(ifs, line);
                        stream.str(line);
                        stream.clear();

                        unsigned int column_index = 0;
                        while (!stream.eof() && (column_index < number_of_columns))
                        {
                            stream >> objectives[objective_index](row_index, column_index);
                            ++column_index;
                        }

                        if (column_index != number_of_columns)
                        {
                            throw std::runtime_error("Not enough data.");
                        }
                    }

                    ++objective_index;
                }
            }
                    
            ifs.close();
        }
    }
}
