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

//EMTGv9 ParallelShootingPhase base class
//Jacob Englander 2-21-2018

#pragma once

#include "phase.h"

#include "PropagatorBase.h"

#include "ParallelShootingStep.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG
{
    namespace Phases
    {
        class ParallelShootingStep;

        class ParallelShootingPhase : public phase
        {
        public:
            //constructor
            ParallelShootingPhase();
            ParallelShootingPhase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            virtual ~ParallelShootingPhase();

            //clone
            virtual ParallelShootingPhase* clone() const = 0;

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            void output_STMs();

            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            //calcbounds goes in the specialized phase
            void setup_calcbounds(std::vector<double>* Xupperbounds,
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

            virtual void calcbounds();

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);


            inline math::Matrix<doubleType>& get_StateAfterInitialCoast() { return this->state_after_initial_coast; }
            inline math::Matrix<doubleType>& get_StateBeforeTerminalCoast() { return this->state_before_terminal_coast; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_state_after_initial_coast() { return this->Derivatives_of_state_after_initial_coast; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_state_before_terminal_coast() { return this->Derivatives_of_state_before_terminal_coast; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_state_after_initial_coast_wrt_Time() { return this->Derivatives_of_state_after_initial_coast_wrt_Time; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_state_before_terminal_coast_wrt_Time() { return this->Derivatives_of_state_before_terminal_coast_wrt_Time; }

            inline size_t get_num_steps() const { return this->num_steps; }

        protected:
            //calcbounds
            void calcbounds_virtual_propellant_tanks();

            virtual void calcbounds_phase_main();

            virtual void calcbounds_initial_coast();

            virtual void calcbounds_terminal_coast();

            virtual void calcbounds_deltav_contribution() {};

            void calcbounds_fixed_inertial_control();

            //process
            virtual void process_phase_main(const std::vector<doubleType>& X,
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
                const bool& needG);

            virtual void process_phase_initial_coast(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_phase_terminal_coast(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_fixed_inertial_control(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            size_t num_steps;
            Astrodynamics::PropagatorBase* InitialCoastPropagatorObject;
            Astrodynamics::PropagatorBase* TerminalCoastPropagatorObject;
            double ForcedCoast_dStepSize_dPropagationVariable;
            double dStepTime_dPhaseFlightTime;

            size_t total_number_of_states_to_integrate;

            boost::ptr_vector<ParallelShootingStep> mySteps;

            //forced coasts
            math::Matrix<doubleType> state_after_initial_coast;
            math::Matrix<double> STM_initial_coast;
            math::Matrix<double> STM_Augmented_initial_coast;
            math::Matrix<double> InitialCoast_dStatedIndependentVariable;
            math::Matrix<doubleType> state_before_terminal_coast;
            math::Matrix<double> STM_terminal_coast;
            math::Matrix<double> STM_Augmented_terminal_coast;
            math::Matrix<double> TerminalCoast_dStatedIndependentVariable;

            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_after_initial_coast;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_before_terminal_coast;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_after_initial_coast_wrt_Time;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_before_terminal_coast_wrt_Time;//Xindex, stateIndex, derivative value

            //derivatives of boundaries
            std::vector<size_t> ListOfVariablesAffectingLeftBoundary;
            std::vector<size_t> ListOfVariablesAffectingRightBoundary;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByVariable; //localXindex, dIndex
            std::vector<size_t> ListOfTimeVariablesAffectingLeftBoundary;
            std::vector<size_t> ListOfTimeVariablesAffectingRightBoundary;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByTimeVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByTimeVariable; //localXindex, dIndex

            //derivatives of state after/before forced coasts
            std::vector< std::vector<size_t> > DerivativesOfStateAfterInitialCoastByVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfStateBeforeTerminalCoastByVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfStateAfterInitialCoastByTimeVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfStateBeforeTerminalCoastByTimeVariable; //localXindex, dIndex

            //tanks
            size_t Xindex_virtual_chemical_fuel_tank;
            size_t Xindex_virtual_electric_propellant_tank;

            //fixed inertial control
            std::vector< std::vector<size_t> > Gindices_fixed_inertial_control_wrt_this_step_control;//step, control variable (num_steps - 1 long)
            std::vector< std::vector<size_t> > Gindices_fixed_inertial_control_wrt_first_step_control;//step, control variable (num_steps - 1 long)
        };
    }//close namespace Phases
}//close namespace EMTG