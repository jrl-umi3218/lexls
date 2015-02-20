#ifndef WORKING_SET
#define WORKING_SET

#include <typedefs.h>

namespace LexLS
{
    namespace internal
    {    
        /** 
            \brief Definition of a working set
        */
        class WorkingSet
        {
        public:
        
            /** 
                \brief Resize and initialize the working set
            
                \param[in] dim Number of constraints
            */
            void resize(Index dim)
            {
                // initialize all constraints as inactive
                all_type.resize(dim,CTR_INACTIVE);
            
                inactive.resize(dim); 
                for (Index CtrIndex=0; CtrIndex<dim; CtrIndex++)
                    inactive[CtrIndex] = CtrIndex;
            }
       
            /**
               \brief Includes in the working set the constraint with index CtrIndex (and sets its type)

               \param[in] CtrIndex Index of constraint in a given LexLSI objective
               \param[in] type     Type of the constraint to be included in the working set

               \note The order of indexes of remaining inactive constraints is modified (but this is not important)
            */                                        
            void activate(Index CtrIndex, ConstraintActivationType type)
            {
                if (all_type[CtrIndex] != CTR_INACTIVE)
                    throw Exception("Cannot activate an active constraint");
            
                std::vector<Index>::iterator it;
                Index ind;
            
                // ----------------------------------------------------------------------------
                // remove the constraint with index CtrIndex from the set of inactive constraints
                // ----------------------------------------------------------------------------
                it = std::find(inactive.begin(), inactive.end(), CtrIndex);
                ind = std::distance(inactive.begin(), it); // it - inactive.begin()
            
                inactive[ind] = inactive.back();
                inactive.pop_back();
            
                // ----------------------------------------------------------------------------
                // add to set of active constraints
                // ----------------------------------------------------------------------------
                all_type[CtrIndex] = type;
                active.push_back(CtrIndex);

                active_ctr_type.push_back(type); // set the type
            }
    
            /** 
                \brief Removes from the working set the constraint with index CtrIndexActive

                \param[in] CtrIndexActive Index of constraint in the working set in a given objective,
                i.e., Obj[ObjIndex].WorkingSet.active[CtrIndexActive] will be removed.
            
                \note The order of indexes of the remaining active constraints is preserved (this would
                be important when we start doing updates).

                \note Note the difference with the function #activate(...).
            */
            void deactivate(Index CtrIndexActive)
            {
                // CtrIndex is the index of WorkingSet.active[CtrIndexActive] in the LexLSI objective
                Index CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex

                if (all_type[CtrIndex] == CTR_INACTIVE)
                    throw Exception("Cannot deactivate an inactive constraint");

                // ----------------------------------------------------------------------------
                // remove the constraint with index CtrIndexActive in the working set
                // ----------------------------------------------------------------------------
                active.erase(active.begin() + CtrIndexActive); // we want to preserve the order of the remaining constraints
                active_ctr_type.erase(active_ctr_type.begin() + CtrIndexActive);

                all_type[CtrIndex] = CTR_INACTIVE;
                // ----------------------------------------------------------------------------
                // add to set of inactive constraints
                // ----------------------------------------------------------------------------
                inactive.push_back(CtrIndex);
            }

            /**
               \brief Returns the number of active constraints
            */
            Index getActiveCtrCount() const
            {
                return active.size();
            }

            /**
               \brief Returns the index of the k-th active constraint
            */                                        
            Index getActiveCtrIndex(Index k) const
            {
                return active[k];
            }

            /**
               \brief Returns the type of the k-th active constraint
            */
            ConstraintActivationType getActiveCtrType(Index k) const
            {
                return active_ctr_type[k];
            }

            /**
               \brief Returns the type of the k-th constraint
            */
            ConstraintActivationType getCtrType(Index k) const
            {
                return all_type[k];
            }

            /**
               \brief Returns the number of inactive constraints
            */
            Index getInactiveCtrCount() const
            {
                return inactive.size();
            }

            /**
               \brief Returns the index of the k-th inactive constraint
            */                                        
            Index getInactiveCtrIndex(Index k) const
            {
                return inactive[k];
            }

            /**
               \brief Returns true if the k-th constraint is active, otherwise returns false
            */                                        
            bool isActive(Index k) const
            {
                if (all_type[k] == CTR_INACTIVE)
                    return false;
                else
                    return true;
            }

            /**
               \brief Prints the contents of the working set
            */                                        
            void print() const
            {
                // -----------------------------------------------------------
                // all
                // -----------------------------------------------------------
                std::cout << " all_type = {";
                std::copy(all_type.begin(),
                          all_type.end(),
                          std::ostream_iterator<ConstraintActivationType>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------
                // type
                // -----------------------------------------------------------
                std::cout << "    type = {";
                std::copy(active_ctr_type.begin(), 
                          active_ctr_type.end(),
                          std::ostream_iterator<ConstraintActivationType>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------
                // active
                // -----------------------------------------------------------
                std::cout << "  active = {";
                std::copy(active.begin(), 
                          active.end(),
                          std::ostream_iterator<Index>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------                
                // inactive
                // -----------------------------------------------------------
                std::cout << "inactive = {";
                std::copy(inactive.begin(),
                          inactive.end(),
                          std::ostream_iterator<Index>(std::cout, " "));
                std::cout << "}" << std::endl;               

                std::cout << std::endl;               
            }

        private:

            /**
               \brief Indexes of active constraints
            */                                        
            std::vector<Index> active;

            /**
               \brief Indexes of inactive constraints
            */                                        
            std::vector<Index> inactive;

            /**
               \brief Type of the active constraints
            */                                        
            std::vector<ConstraintActivationType> active_ctr_type;

            /**
               \brief Specifies the type of all constraints

               \note Introduced for convenience
            */                                        
            std::vector<ConstraintActivationType> all_type;
        };

    } // END namespace internal
    
} // END namespace LexLS

#endif // WORKING_SET
