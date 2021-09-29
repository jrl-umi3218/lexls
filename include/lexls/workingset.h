/*
 * Copyright 2013-2021 INRIA
 */

#pragma once

#include <lexls/typedefs.h>

#include <numeric> // for std::iota

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
            inline void resize(Index dim)
            {
                all_type.resize(dim);
                active.reserve(dim);
                active_ctr_type.reserve(dim);
                inactive.resize(dim);
                reset();
            }

            /**
                \brief Reset the working set to all constraints being inactive.
            */
            inline void reset() {
                // initialize all constraints as inactive
                std::fill(all_type.begin(), all_type.end(), CTR_INACTIVE);
                active.clear();
                active_ctr_type.clear();
                // initialize inactive with [0, 1, 2, ...]
                std::iota(inactive.begin(), inactive.end(), 0);
            }

            /**
               \brief Includes in the working set the constraint with index CtrIndex (and sets its type)

               \param[in] CtrIndex Index of constraint in a given LexLSI objective
               \param[in] type     Type of the constraint to be included in the working set

               \note The order of indexes of remaining inactive constraints is modified (but this is not important)
            */
            inline void activate(Index CtrIndex, ConstraintActivationType type)
            {
                if (all_type[CtrIndex] != CTR_INACTIVE)
                {
                    throw Exception("Cannot activate an active constraint");
                }

                // ----------------------------------------------------------------------------
                // remove the constraint with index CtrIndex from the set of inactive constraints
                // ----------------------------------------------------------------------------
                Index ind = getCtrIndex(CtrIndex); // get the index among the inactive constraints

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
            inline void deactivate(Index CtrIndexActive)
            {
                // CtrIndex is the index of WorkingSet.active[CtrIndexActive] in the LexLSI objective
                Index CtrIndex = getActiveCtrIndex(CtrIndexActive); // CtrIndexActive --> CtrIndex

                // \todo This test is useless: attempting to deactivate a non activated constraint will
                // make the above line crash.
                if (all_type[CtrIndex] == CTR_INACTIVE)
                {
                    throw Exception("Cannot deactivate an inactive constraint");
                }

                // ----------------------------------------------------------------------------
                // remove the constraint with index CtrIndexActive in the working set
                // ----------------------------------------------------------------------------
                active.erase(active.begin()
                             + CtrIndexActive); // we want to preserve the order of the remaining constraints
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
            inline Index getActiveCtrCount() const
            {
                return static_cast<Index>(active.size());
            }

            /**
               \brief Returns the index of the k-th active constraint
            */
            inline Index getActiveCtrIndex(Index k) const
            {
                return active[k];
            }

            /**
               \brief Returns the type of the k-th active constraint
            */
            inline ConstraintActivationType getActiveCtrType(Index k) const
            {
                return active_ctr_type[k];
            }

            /**
               \brief Returns the type of the k-th constraint
            */
            inline ConstraintActivationType getCtrType(Index k) const
            {
                return all_type[k];
            }

            /**
               \brief If the constraint with index k is active then return its index among the
               active constraints. If the constraint with index k is inactive then return its index
               among the inactive constraints.
            */
            inline Index getCtrIndex(Index k) const
            {
                if (isActive(k))
                {
                    auto it = std::find(active.begin(), active.end(), k);

                    return static_cast<Index>(std::distance(active.begin(), it)); // it - active.begin()
                }
                else
                {
                    auto it = std::find(inactive.begin(), inactive.end(), k);

                    return static_cast<Index>(std::distance(inactive.begin(), it)); // it - inactive.begin()
                }
            }

            /**
               \brief Returns the number of inactive constraints
            */
            inline Index getInactiveCtrCount() const
            {
                return static_cast<Index>(inactive.size());
            }

            /**
               \brief Returns the index of the k-th inactive constraint
            */
            inline Index getInactiveCtrIndex(Index k) const
            {
                return inactive[k];
            }

            /**
               \brief Returns true if the k-th constraint is active, otherwise returns false
            */
            inline bool isActive(Index k) const
            {
                if (all_type[k] == CTR_INACTIVE)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }

            /**
               \brief Prints the contents of the working set
            */
            inline void print() const
            {
                // -----------------------------------------------------------
                // all
                // -----------------------------------------------------------
                std::cout << " all_type = {";
                std::copy(all_type.begin(), all_type.end(),
                          std::ostream_iterator<ConstraintActivationType>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------
                // type
                // -----------------------------------------------------------
                std::cout << "     type = {";
                std::copy(active_ctr_type.begin(), active_ctr_type.end(),
                          std::ostream_iterator<ConstraintActivationType>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------
                // active
                // -----------------------------------------------------------
                std::cout << "   active = {";
                std::copy(active.begin(), active.end(), std::ostream_iterator<Index>(std::cout, " "));
                std::cout << "}" << std::endl;

                // -----------------------------------------------------------
                // inactive
                // -----------------------------------------------------------
                std::cout << " inactive = {";
                std::copy(inactive.begin(), inactive.end(), std::ostream_iterator<Index>(std::cout, " "));
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
