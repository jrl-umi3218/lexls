/*
 * Copyright 2013-2021 INRIA
 */

#pragma once

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
                previous_operation(OPERATION_UNDEFINED)
            {
                previous_ctr_identifier.set(0,0,CTR_INACTIVE);
            }

            TerminationStatus update(OperationType operation,
                                     ConstraintIdentifier ctr_identifier,
                                     std::vector<Objective> &Obj,
                                     bool &cycling_detected)
            {
                cycling_detected = false;
                if ((operation == OPERATION_ADD) &&
                    (previous_operation == OPERATION_REMOVE))
                {
                    if (ctr_identifier == previous_ctr_identifier)
                    {
                        if (counter >= max_counter)
                        {
                            return PROBLEM_SOLVED_CYCLING_HANDLING;
                        }
                        else
                        {
                            relax_bounds(Obj);
                            cycling_detected = true;
                        }
                    }
                }
                previous_operation      = operation;
                previous_ctr_identifier = ctr_identifier;

                return TERMINATION_STATUS_UNKNOWN;
            }

            void relax_bounds(std::vector<Objective> &Obj)
            {
                Obj[previous_ctr_identifier.obj_index]
                    .relax_bounds(previous_ctr_identifier.ctr_index,
                                  previous_ctr_identifier.ctr_type,
                                  relax_step);

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
               \brief Constraint involved in the previous operation
            */
            ConstraintIdentifier previous_ctr_identifier;
        };

    } // END namespace internal

} // END namespace LexLS
