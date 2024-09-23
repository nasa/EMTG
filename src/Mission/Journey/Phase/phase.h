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

//EMTGv9 phase base class

//#ifndef SIZE_MAX
//#define SIZE_MAX 65535
//#endif

#pragma once

#include "doubleType.h"

#include <string>
#include <vector>
#include <fstream>

#include "missionoptions.h"
#include "universe.h"
#include "LaunchVehicle.h"
#include "Spacecraft.h"
#include "writey_thing.h"
#include "sparsey_thing.h"

#include "arrival.h"
#include "departure.h"

namespace EMTG
{
    //forward declare arrival and departure
    class DepartureEvent;
    class ArrivalEvent;

    namespace Phases
    {
        class phase : public writey_thing, public sparsey_thing
        {
        public:
            //constructor
            phase();

            phase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions,
                const size_t& numStatesToPropagate,
                const size_t& numMatchConstraints);

            //clone
            virtual phase* clone() const = 0;

            //destructor
            virtual ~phase();

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount) = 0; 
            
            void output_initial_TCM(std::ofstream& outputfile,
                size_t& eventcount);
            
            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream & acceleration_model_file) = 0;

            virtual void output_STMs() = 0;

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget) = 0;


            //set
            void setName(const std::string& name) { this->name = name; }

            //get
            inline size_t getX_index_of_first_decision_variable_in_this_phase() const { return this->X_index_of_first_decision_variable_in_this_phase; }
            inline double get_initial_TCM_magnitude() const { return this->initial_TCM_magnitude; }
            inline BoundaryEvents::DepartureEvent* getDepartureEvent() const { return this->myDepartureEvent; }
            inline BoundaryEvents::ArrivalEvent* getArrivalEvent() const { return this->myArrivalEvent; }
            inline double getsynodic_period() const { return this->synodic_period; }
            inline double getPhaseDutyCycle() const { return this->PhaseDutyCycle; }
            PropagatorType getPropagatorType() const { return this->myPropagatorType; }

            doubleType getDeterministicDeltav() const { return this->PhaseTotalDeterministicDeltav; }
            doubleType getStatisticalDeltav() const { return this->PhaseTotalStatisticalDeltav; }
            double getInitialCoastDuration() const { return this->InitialCoastDuration; }
            double getTerminalCoastDuration() const { return this->TerminalCoastDuration; }
            doubleType getFinalMass() const { return this->myArrivalEvent->get_state_after_event()(6); }
            bool getIsLastPhaseInJourney() const { return this->isLastPhaseInJourney; }

            inline std::vector<std::string>& get_matchPointConstraintNames() { return this->matchPointConstraintNames; }
            inline std::vector<size_t>& get_matchPointConstraintStateIndex() { return this->matchPointConstraintStateIndex; }
            inline math::Matrix<double>& get_continuity_constraint_scale_factors() { return this->continuity_constraint_scale_factors; }
            inline size_t get_numMatchConstraints() { return this->numMatchConstraints; }

            //calcbounds
            virtual void setup_calcbounds(std::vector<double>* Xupperbounds,
                std::vector<double>* Xlowerbounds,
                std::vector<double>* X_scale_factors,
                std::vector<double>* Fupperbounds,
                std::vector<double>* Flowerbounds,
				std::vector<double>* F_scale_factors,
                std::vector<std::string>* Xdescriptions,
                std::vector<std::string>* Fdescriptions,
                std::vector<size_t>* iGfun,
                std::vector<size_t>* jGvar,
                std::vector<std::string>* Gdescriptions,
                std::vector<size_t>* iAfun,
                std::vector<size_t>* jAvar,
                std::vector<std::string>* Adescriptions,
                std::vector<double>* A);

            //calcbounds goes in the specialized phase
            virtual void calcbounds() = 0;

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

        protected:
            void reset_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void calcbounds_phase_left_boundary();

            virtual void calcbounds_phase_flight_time();

            void calcbounds_phase_right_boundary();

            virtual void calcbounds_phase_main() = 0;

            virtual void calcbounds_virtual_propellant_tanks() = 0;

            virtual void calcbounds_deltav_contribution() = 0;

            virtual void process_phase_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            void process_phase_left_boundary(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_phase_flight_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_phase_right_boundary(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG); 

            virtual void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;


            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            HardwareModels::Spacecraft* mySpacecraft;
            HardwareModels::LaunchVehicle* myLaunchVehicle;
            JourneyOptions* myJourneyOptions;
            BoundaryEvents::ArrivalEvent* PreviousPhaseArrivalEvent;
            size_t X_index_of_first_decision_variable_in_this_phase;
            size_t F_index_of_first_constraint_in_this_phase;
            std::vector<std::string> stateVectorNames;
            std::vector<std::string> matchPointConstraintNames;
            std::vector<size_t> matchPointConstraintStateIndex;
            bool isFirstPhaseInMission;
            bool isFirstPhaseInJourney;
            bool isLastPhaseInJourney;
            PropagatorType myPropagatorType;

            size_t numMatchConstraints;
            size_t numStatesToPropagate;
            math::Matrix<double> continuity_constraint_scale_factors;

            //pointers to boundary events
            BoundaryEvents::DepartureEvent* myDepartureEvent;
            BoundaryEvents::ArrivalEvent* myArrivalEvent;

            //boundary states
            math::Matrix<doubleType> state_at_beginning_of_phase;
            math::Matrix<doubleType> state_at_end_of_phase;

            //truth table relating encoded states to match constraints - if an entry is zero, that derivative is zero and need not be computed
            std::vector< std::vector<size_t> > TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates;
            std::vector< std::vector<size_t> > TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates;
            std::vector< size_t > TruthTable_MatchConstraints_Derivative_wrt_PropagationTime;
            size_t stateIndex_phase_propagation_variable;

            //*********************delta-v and derivatives
            size_t Findex_VirtualPhaseTotalDeterministicDeltav;
            doubleType VirtualPhaseTotalDeterministicDeltav;
            doubleType PhaseTotalDeterministicDeltav;
            doubleType PhaseTotalStatisticalDeltav;

            //staging
            bool StageBeforePhase; //if true, then we need to apply a staging mass constraint
            size_t Gindex_StageMass_LeftBoundaryMass;

            //propellant used
            doubleType virtual_electric_propellant_used;
            doubleType virtual_chemical_fuel_used;
            doubleType virtual_chemical_oxidizer_used;
            doubleType electric_propellant_used;
            doubleType chemical_fuel_used;
            doubleType chemical_oxidizer_used;
            bool hasElectricManeuver;
            bool hasMonopropManeuver;
            bool hasBipropManeuver;

            size_t Findex_VirtualElectricPropellantConstraint;
            size_t Findex_VirtualChemicalFuelConstraint;
            size_t Findex_VirtualChemicalOxidizerConstraint;

            size_t Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant;
            size_t Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
            size_t Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer;

            //TCM stuff
            bool hasInitialTCM;
            double initial_TCM_magnitude;
            math::Matrix<doubleType> state_after_initial_TCM;
            math::Matrix<double> TCMTM; //TCM transition matrix

            //forced coast stuff
            bool hasInitialCoast;
            double InitialCoastDuration;
            bool hasTerminalCoast;
            double TerminalCoastDuration;

            //time
            doubleType PhaseFlightTime;
            double synodic_period;
            double integration_step_length;
            size_t Xindex_PhaseFlightTime;
            doubleType LaunchDate; //needed for power things
            size_t Xindex_LaunchDate;

            //propulsion things - all low-thrust phases have a duty cycle
            double PhaseDutyCycle;

            //plotting things
            double EphemerisOutputResolution;

            //output
            math::Matrix<doubleType> output_state;
            math::Matrix<doubleType> temp_state;
            doubleType output_power;
            doubleType output_active_power;
        };
        
        inline phase * new_clone(phase const & other)
        {
            return other.clone();
        }
    }//namespace Phases
} /* namespace EMTG */
