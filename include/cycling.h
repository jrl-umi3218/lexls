// Time-stamp: <2014-09-17 11:46:10 drdv>
#ifndef CYCLING
#define CYCLING

/**
   \brief Cycling detection
*/
namespace LexLS
{   
    /** 
        \brief A class used for cycling detection.
    */
    class CyclingHandlerType
    {    
    public:
    
        CyclingHandlerType():
            counter(0),
            max_counter(50),
            relax_step(1e-08),
            previous_operation(UNDEFINED)
        {

        }

        void print_counter()
        {
            printf("   counter = %d\n",counter);
        }

        void print()
        {
            printf("============================================= \n");
            printf("CtrAdded (%d): \n", (Index)CtrAdded.size());
            for (Index k=0; k<CtrAdded.size(); k++)
            {
                printf("     ");
                CtrAdded[k].print();
            }

            printf("============================================= \n");
            printf("CtrRemoved (%d): \n", (Index)CtrRemoved.size());
            for (Index k=0; k<CtrRemoved.size(); k++)
            {
                printf("     ");
                CtrRemoved[k].print();
            }
        }

        TerminationStatus update(OperationType operation, 
                                 ConstraintIdentifierType ctr_identifier, 
                                 std::vector<Objective> &Obj,
                                 Index iter,
                                 bool print_flag=false)
        {
            // this is executed only the first time this function is called
            if (previous_operation == UNDEFINED)
            {
                if (operation == ADD_CONSTRAINT)
                    previous_operation = REMOVE_CONSTRAINT;
                else if (operation == REMOVE_CONSTRAINT)
                    previous_operation = ADD_CONSTRAINT;
            }

            if (print_flag)
            {
                printf("============================================= \n");
                printf("iter = %d, counter = %d\n   previous_operation = %s\n   operation          = %s \n", 
                       iter,
                       counter,
                       previous_operation == ADD_CONSTRAINT ? "ADD" : "REMOVE", 
                       operation == ADD_CONSTRAINT ? "ADD" : "REMOVE");
            }

            if (operation == ADD_CONSTRAINT)
            {
                if (previous_operation == REMOVE_CONSTRAINT)
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

                    previous_operation = ADD_CONSTRAINT;
                    CtrAdded.clear();
                }
                
                CtrAdded.push_back(ctr_identifier);
            }
            else if (operation == REMOVE_CONSTRAINT)
            {
                if (previous_operation == ADD_CONSTRAINT)
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
   
                    previous_operation = REMOVE_CONSTRAINT;
                    CtrRemoved.clear();
                }

                CtrRemoved.push_back(ctr_identifier);
            }    

            if (print_flag)
                print();

            return TERMINATION_STATUS_UNKNOWN;
        }

        bool check_condition()
        {
            if (CtrRemoved.size() == 0 || CtrAdded.size() == 0) // disregard cases with empty containers 
                return false;

            bool condition = true;
            if (CtrRemoved.size() == CtrAdded.size())
            {
                for (Index k=0; k<CtrRemoved.size(); k++)
                {
                    if (std::find(CtrRemoved.begin(), CtrRemoved.end(), CtrAdded[k]) == CtrRemoved.end())
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
            for (Index k=0; k<CtrRemoved.size(); k++)
            {
                Obj[CtrRemoved[k].getObjIndex()]
                    .relax_bounds(CtrRemoved[k].getCtrIndex(), CtrRemoved[k].getCtrType(), relax_step);
            }
            counter++;
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
        std::vector<ConstraintIdentifierType> CtrAdded;
        
        /** 
            \brief List of removed constraints after the last added constraints
        */       
        std::vector<ConstraintIdentifierType> CtrRemoved;        
    };
    
} // END namespace LexLS

#endif // CYCLING

