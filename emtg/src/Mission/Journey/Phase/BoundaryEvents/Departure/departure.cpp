// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

#include "departure.h"

#include "EMTG_enums.h"
#include "orbit_element_conversions.h"
#include "KeplerPropagatorTimeDomain.h"
#include "SpecializedBoundaryConstraintFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        DepartureEvent::DepartureEvent(const std::string& name,
                                       const size_t& journeyIndex,
                                       const size_t& phaseIndex,
                                       size_t& stageIndex,
                                       Astrodynamics::universe* Universe,
                                       HardwareModels::Spacecraft* mySpacecraft,
                                       missionoptions* myOptions,
                                       ArrivalEvent* PreviousPhaseArrivalEvent)
            
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions,
                             PreviousPhaseArrivalEvent);
        }//end constructor
        
        
        
        void DepartureEvent::initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->PreviousPhaseArrivalEvent = PreviousPhaseArrivalEvent;
            this->hasFixedInitialMass = false;
            this->initial_mass_increment = 0.0;

            if (this->phaseIndex == 0)
            {
                this->isFirstEventInJourney = true;

                if (this->journeyIndex == 0)
                {
                    this->isFirstEventInMission = true;
                    this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                }
            }
            else
                this->isFirstEventInJourney = false;

            //boundary constraints
            this->construct_boundary_constraints();
        }//end initialize()
        
        void DepartureEvent::construct_boundary_constraints(std::vector<std::string> givenConstraints)
        {
            //first construct this event's tag
            std::string Tag = "p" + std::to_string(this->phaseIndex)
                + "_departure";

            std::string TagEnd = "pEnd_departure";

            //clear the current constraint vector
            this->mySpecializedConstraints.clear();

            //are we adding constraints from journeyOptions or from an input vector?
            std::vector<std::string>* constraintsToAdd;

            if (!givenConstraints.empty())
            {
                constraintsToAdd = &givenConstraints;
            }
            else
            {
                constraintsToAdd = &this->myJourneyOptions->BoundaryConstraintDefinitions;
            }

            //now, loop over constraints to see if they are relevant
            for (std::string& constraint : *constraintsToAdd)
            {
                if (constraint.find("#") != 0) //don't create a constraint if it is commented out
                {
                    if (constraint.find(Tag) < 1024 || (constraint.find(TagEnd) < 1024 && this->phaseIndex == this->myJourneyOptions->number_of_phases - 1))
                    {
                        if (constraint.find("monoprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Monoprop;
                        }
                        else if (constraint.find("biprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Biprop;
                        }
                        else
                        {
                            this->mySpecializedConstraints.push_back(BoundaryEvents::SpecializedConstraints::create_boundary_event_constraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                (BoundaryEventBase*)this,
                                this->myUniverse,
                                this->mySpacecraft,
                                this->myOptions,
                                constraint,
                                "departure"));
                        }
                    }
                }//end constraint comment check
            }//end loop over constraints
        }//end construct_boundary_constraints
        
        //******************************************calcbounds methods
        void DepartureEvent::calcbounds_event_left_side()
        {
            //wait time and mass continuity
            if (this->hasWaitTime)
            {
                Xlowerbounds->push_back(this->myJourneyOptions->wait_time_bounds[0]);
                Xupperbounds->push_back(this->myJourneyOptions->wait_time_bounds[1]);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(prefix + "wait time");
                this->Xindices_EventLeftEpoch.push_back(this->Xdescriptions->size() - 1);
            }
        }//end calcbounds_event_left_side()

        void DepartureEvent::calcbounds_mass_multipliers()
        {
            if (this->phaseIndex == 0
                && this->myJourneyOptions->variable_mass_increment)
            {
                Xlowerbounds->push_back(math::SMALL);
                Xupperbounds->push_back(1.0);
                X_scale_factors->push_back(1.0);
                Xdescriptions->push_back(this->prefix + "journey initial mass increment multiplier");

                this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple<size_t, size_t, double>(Xdescriptions->size() - 1, 6, 1.0));
                this->derivative_index_of_journey_initial_mass_increment_multiplier = this->Derivatives_of_StateBeforeEvent.size() - 1;
            }

            //journey initial mass constraint
            if (this->myJourneyOptions->constrain_initial_mass)
            {
                this->Flowerbounds->push_back(-this->myJourneyOptions->maximum_initial_mass / this->myJourneyOptions->maximum_mass);
                this->Fupperbounds->push_back(0.0);
                this->Fdescriptions->push_back(this->prefix + "journey initial mass constraint");

                //has dependencies on all variables that affect initial mass
                for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforeEvent[dIndex]);

                    if (stateIndex == 6)
                    {
                        this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables.push_back(dIndex);

                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforeEvent[dIndex]);

                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables);
                    }
                }

                //time variables
                for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                    if (stateIndex == 6)
                    {
                        this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time.push_back(dIndex);

                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time);
                    }
                }
            }//end journey initial mass constraint
        }//end calcbounds_mass_multipliers()

        void DepartureEvent::calcbounds_left_mass_continuity_constraint()
        {
            //Step 1: make the constraint
            Flowerbounds->push_back(-1.0e-13);
            Fupperbounds->push_back(1.0e-13);
            Fdescriptions->push_back(prefix + "left mass continuity constraint");

            //derivatives with respect to previous phase state after arrival event
            std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent
                = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();

            std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time
                = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //Step 2: we need to store the dIndex of, and create a sparsity entry for, anything that affects the previous event's "mass after event"
            for (size_t dIndex = 0; dIndex < PreviousEvent_Derivatives_of_StateAfterEvent.size(); ++dIndex)
            {
                size_t stateIndex = std::get<1>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex]);
if (stateIndex == 6)
{
    this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables.push_back(dIndex);
    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
        std::get<0>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex]),
        this->Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables);
}
            }//end loop over decision variables that affect the previous event's "state after event"

            //with respect to time
            for (size_t dIndex = 0; dIndex < PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
            {
                size_t stateIndex = std::get<1>(PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);
                if (stateIndex == 6)
                {
                    this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        std::get<0>(PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                        this->Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time);
                }
            }//end loop over decision variables that affect the previous event's "state after event"


            //Step 3: we need to store the dIndex of, and create a sparsity entry for, anything that affects the current event's "mass before event"
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
            {
                size_t stateIndex = std::get<1>(this->Derivatives_of_StateBeforeEvent[dIndex]);
                if (stateIndex == 6)
                {
                    this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex]),
                        this->Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables);
                }
            }//end loop over decision variables that affect the current event's "state before event"

            //with respect to time
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
            {
                size_t stateIndex = std::get<1>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);
                if (stateIndex == 6)
                {
                    this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        std::get<0>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]),
                        this->Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time);
                }
            }//end loop over decision variables that affect the previous event's "state after event"
        }//end calcbounds_left_mass_continuity_constraint()

        //**************************************process functions
        void DepartureEvent::
            process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            if (this->hasWaitTime)
                this->EventWaitTime = X[Xindex++];
        }//end process_event_left_side()


        void DepartureEvent::
            process_mass_multipliers(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            if (this->isFirstEventInJourney
                && this->myJourneyOptions->variable_mass_increment)
            {
                this->journey_initial_mass_increment_multiplier = X[Xindex++];
                double mass_scale = (this->myJourneyOptions->maximum_starting_mass_increment
                    - this->myJourneyOptions->minimum_starting_mass_increment);

                this->initial_mass_increment = this->journey_initial_mass_increment_multiplier
                    * mass_scale + this->myJourneyOptions->minimum_starting_mass_increment;


                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->derivative_index_of_journey_initial_mass_increment_multiplier])
                    = mass_scale;
            }
            else
            {
                this->initial_mass_increment = this->myJourneyOptions->fixed_starting_mass_increment;
            }

            this->state_before_event(6) += this->initial_mass_increment;

            //journey initial mass constraint
            if (this->myJourneyOptions->constrain_initial_mass)
            {
                F[Findex++] = (this->state_before_event(6) - this->myJourneyOptions->maximum_initial_mass) / this->myJourneyOptions->maximum_mass;

                //has dependencies on all variables that affect initial mass
                for (size_t varIndex = 0; varIndex < this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables.size(); ++varIndex)
                {
                    size_t dIndex = this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables[varIndex];

                    double TheDerivative = std::get<2>(Derivatives_of_StateBeforeEvent[dIndex]);

                    size_t Gindex = this->Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * TheDerivative
                        / this->myJourneyOptions->maximum_mass;
                }//end non-time derivatives

                //time derivatives
                for (size_t varIndex = 0; varIndex < this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time.size(); ++varIndex)
                {
                    size_t dIndex = this->dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time[varIndex];

                    double TheDerivative = std::get<2>(Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                    size_t Gindex = this->Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * TheDerivative
                        / this->myJourneyOptions->maximum_mass;
                }//end time derivatives
            }//end journey initial mass constraint
        }//end process_mass_multipliers()



        void DepartureEvent::process_left_mass_continuity_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: constraint
            math::Matrix<doubleType>& PreviousEventStateAfterEvent = this->PreviousPhaseArrivalEvent->get_state_after_event();

            F[Findex++] = (this->state_before_event(6) - PreviousEventStateAfterEvent(6)) * this->myUniverse->continuity_constraint_scale_factors(6);

            //derivatives
            if (needG)
            {
                //derivatives with respect to previous phase state after arrival event
                std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent
                    = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();

                std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time
                    = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //Step 2: we need to store the dIndex of, and create a sparsity entry for, anything that affects the previous event's "mass after event"
                for (size_t varIndex = 0; varIndex < this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables.size(); ++varIndex)
                {
                    size_t Gindex = this->Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    size_t dIndex = this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables[varIndex];
                    double Gentry = std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex]);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }//end loop over decision variables that affect the previous event's "state after event"

                //with respect to time
                for (size_t varIndex = 0; varIndex < this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time.size(); ++varIndex)
                {
                    size_t Gindex = this->Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    size_t dIndex = this->dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time[varIndex];
                    double Gentry = std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }//end loop over decision variables that affect the previous event's "state after event"


                //Step 3: we need to store the dIndex of, and create a sparsity entry for, anything that affects the current event's "mass before event"
                for (size_t varIndex = 0; varIndex < this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables.size(); ++varIndex)
                {
                    size_t Gindex = this->Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    size_t dIndex = this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables[varIndex];
                    double Gentry = std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }//end loop over decision variables that affect the previous event's "state after event"

                 //with respect to time
                for (size_t varIndex = 0; varIndex < this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time.size(); ++varIndex)
                {
                    size_t Gindex = this->Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    size_t dIndex = this->dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time[varIndex];
                    double Gentry = std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }//end loop over decision variables that affect the previous event's "state after event"
            }//end derivatives
        }//end process_left_mass_continuity_constraint()

         //******************************************output methods
        void DepartureEvent::output_mass_increment(std::ofstream& outputfile)
        {
            outputfile << "Journey initial mass increment: " << this->initial_mass_increment _GETVALUE << " kg" << std::endl;
        }
        
        void DepartureEvent::output_ephemeris(std::ofstream& outputfile)
        {
            if (this->isFirstEventInMission || this->hasWaitTime)
                this->BoundaryEventBase::output_ephemeris(outputfile);
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG