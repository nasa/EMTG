
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

#pragma once

#include "BoundaryEventBase.h"
#include "arrival.h"
#include "PropagatorBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class DepartureEvent : virtual public BoundaryEventBase
        {
        public:
            //default constructor
            DepartureEvent() : BoundaryEventBase::BoundaryEventBase() {};

            //specialized constructor
            DepartureEvent(const std::string& name,
                           const size_t& journeyIndex,
                           const size_t& phaseIndex,
                           size_t& stageIndex,
                           Astrodynamics::universe* Universe,
                           HardwareModels::Spacecraft* mySpacecraft,
                           missionoptions* myOptions,
                           ArrivalEvent* PreviousPhaseArrivalEvent);

            virtual void initialize(const std::string& name,
                                    const size_t& journeyIndex,
                                    const size_t& phaseIndex,
                                    size_t& stageIndex,
                                    Astrodynamics::universe* Universe,
                                    HardwareModels::Spacecraft* mySpacecraft,
                                    missionoptions* myOptions,
                                    ArrivalEvent* PreviousPhaseArrivalEvent);

            virtual void construct_boundary_constraints(std::vector<std::string> givenConstraints = {});

            //destructor
            virtual ~DepartureEvent() {};

            //output
            
            virtual void output_mass_increment(std::ofstream& outputfile);

            virtual math::Matrix<doubleType> get_periapse_state() = 0;
            virtual void output_periapse_state(size_t& flybyIndex, std::ofstream& outputfile) {};

            inline doubleType getInitialMassIncrement() const { return this->initial_mass_increment; }
            inline bool getHasWaitTime() const { return this->hasWaitTime; }

            virtual void output_ephemeris(std::ofstream& outputfile);

            //this is a stub, but some events may override it with an actual maneuver
            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget) {};

        protected:

            //method to calculate event left side
            virtual void calcbounds_event_left_side(std::vector<size_t> timeVariables) = 0;

            virtual void calcbounds_event_left_side(); //does the wait time

            void calcbounds_mass_multipliers();

            void calcbounds_left_mass_continuity_constraint();

            virtual void calcbounds_event_right_side() = 0;

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            void process_left_mass_continuity_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_mass_multipliers(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_event_right_side(const std::vector<doubleType>& X, //this is just to handle the derivative with respect to wait time
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;


            //fields
            ArrivalEvent* PreviousPhaseArrivalEvent;//pointer to previous phase arrival event
            bool isFirstEventInJourney;
            bool hasFixedInitialMass;

            //times
            doubleType EventWaitTime;
            //I have to initialize this in the header because only the derived class can make it true
            //and the base class's constructor is only ever called AFTER the derived class does this
            bool hasWaitTime = false;

            //mass multipliers
            doubleType journey_initial_mass_increment_multiplier;
            doubleType initial_mass_increment;
            size_t derivative_index_of_journey_initial_mass_increment_multiplier;

            //left mass continuity
            std::vector<size_t> dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables;
            std::vector<size_t> dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables;
            std::vector<size_t> dIndex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time;
            std::vector<size_t> dIndex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time;
            std::vector<size_t> Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables;
            std::vector<size_t> Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables;
            std::vector<size_t> Gindex_LeftMassContinuity_PreviousArrivalEventDecisionVariables_wrt_Time;
            std::vector<size_t> Gindex_LeftMassContinuity_StateBeforeEventDecisionVariables_wrt_Time;

            size_t Xindex_previous_right_mass;
            size_t Xindex_current_left_mass;
            size_t Gindex_left_mass_continuity_wrt_current_event_left_mass;
            size_t Gindex_left_mass_continuity_wrt_previous_event_right_mass;

            //journey initial mass constraint
            std::vector<size_t> dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables;
            std::vector<size_t> dIndex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time;
            std::vector<size_t> Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables;
            std::vector<size_t> Gindex_journey_initial_mass_constraint_wrt_StateBeforeEventDecisionVariables_wrt_Time;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG