#ifndef CYCLING
#define CYCLING

namespace LexLS
{   
    namespace internal
    {
        /** 
            \brief A class used for cycling detection.
        */
        class CyclingHandler
        {    
        public:
    
            CyclingHandler():
                counter(0),
                max_counter(50),
                relax_step(1e-08),
                previous_operation(OPERATION_UNDEFINED){}

            void print_counter()
            {
                printf("   counter = %d\n",counter);
            }

            void print()
            {
                printf("============================================= \n");
                printf("ctr_added (%d): \n", (Index)ctr_added.size());
                for (Index k=0; k<(Index)ctr_added.size(); k++)
                {
                    printf("     ");
                    ctr_added[k].print();
                }

                printf("============================================= \n");
                printf("ctr_removed (%d): \n", (Index)ctr_removed.size());
                for (Index k=0; k<(Index)ctr_removed.size(); k++)
                {
                    printf("     ");
                    ctr_removed[k].print();
                }
            }

            TerminationStatus update(OperationType operation, 
                                     ConstraintIdentifier ctr_identifier, 
                                     std::vector<Objective> &Obj,
                                     Index iter,
                                     bool print_flag=false)
            {
                // this is executed only the first time this function is called
                if (previous_operation == OPERATION_UNDEFINED)
                {
                    if (operation == OPERATION_ADD)
                        previous_operation = OPERATION_REMOVE;
                    else if (operation == OPERATION_REMOVE)
                        previous_operation = OPERATION_ADD;
                }

                if (print_flag)
                {
                    printf("============================================= \n");
                    printf("iter = %d, counter = %d\n   previous_operation = %s\n   operation          = %s \n", 
                           iter,
                           counter,
                           previous_operation == OPERATION_ADD ? "ADD" : "REMOVE", 
                           operation == OPERATION_ADD ? "ADD" : "REMOVE");
                }

                if (operation == OPERATION_ADD)
                {
                    if (previous_operation == OPERATION_REMOVE)
                    {
                        if (check_condition())
                        {
                            if (counter >= max_counter)
                            {
                                if (print_flag)
                                    print();

                                return PROBLEM_SOLVED_CYCLING_HANDLING;
                            }
                            else
                            {
                                relax_bounds(Obj);
                            }
                        }

                        previous_operation = OPERATION_ADD;
                        ctr_added.clear();
                    }
                
                    ctr_added.push_back(ctr_identifier);
                }
                else if (operation == OPERATION_REMOVE)
                {
                    if (previous_operation == OPERATION_ADD)
                    {
                        if (check_condition())
                        {
                            if (counter >= max_counter)
                            {
                                if (print_flag)
                                    print();
                            
                                return PROBLEM_SOLVED_CYCLING_HANDLING;
                            }
                            else
                            {
                                relax_bounds(Obj);
                            }
                        }
   
                        previous_operation = OPERATION_REMOVE;
                        ctr_removed.clear();
                    }

                    ctr_removed.push_back(ctr_identifier);
                }    

                if (print_flag)
                    print();

                return TERMINATION_STATUS_UNKNOWN;
            }

            bool check_condition()
            {
                if (ctr_removed.size() == 0 || ctr_added.size() == 0) // disregard cases with empty containers 
                    return false;

                bool condition = true;
                if (ctr_removed.size() == ctr_added.size())
                {
                    for (Index k=0; k<(Index)ctr_removed.size(); k++)
                    {
                        if (std::find(ctr_removed.begin(), ctr_removed.end(), ctr_added[k]) == ctr_removed.end())
                        {
                            condition = false;
                            break;
                        }
                    }
                }
            
                return condition;
            }

            void relax_bounds(std::vector<Objective> &Obj)
            {
                for (Index k=0; k<(Index)ctr_removed.size(); k++)
                {
                    Obj[ctr_removed[k].getObjIndex()]
                        .relax_bounds(ctr_removed[k].getCtrIndex(), ctr_removed[k].getCtrType(), relax_step);
                }
                counter++;
            }

            void set_max_counter(Index max_counter_)
            {
                max_counter = max_counter_;
            }

            void set_relax_step(RealScalar relax_step_)
            {
                relax_step = relax_step_;
            }

            Index get_counter() const
            {
                return counter;
            }

            /** 
                \brief Number of relaxations performed
            */ 
            Index counter;

            /** 
                \brief Maximum number of relaxations to be performed
            */ 
            Index max_counter;

            /** 
                \brief Relaxation step
            */ 
            RealScalar relax_step;

            /** 
                \brief Operation performed during the previous iteration
            */ 
            OperationType previous_operation;
        
            /** 
                \brief List of added constraints after the last removed constraints
            */       
            std::vector<ConstraintIdentifier> ctr_added;
        
            /** 
                \brief List of removed constraints after the last added constraints
            */       
            std::vector<ConstraintIdentifier> ctr_removed;        
        };

    } // END namespace internal
    
} // END namespace LexLS

#endif // CYCLING

