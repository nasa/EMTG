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

//EMTGv9 TwoPointShootingLowThrustPhase
//Jacob Englander 6-23-2017

#pragma once

#include "TwoPointShootingPhase.h"

namespace EMTG 
{
    namespace Phases
    {
        class TwoPointShootingLowThrustPhase : public TwoPointShootingPhase
        {
        public:
            //constructor
            TwoPointShootingLowThrustPhase() : TwoPointShootingPhase::TwoPointShootingPhase() {};
            TwoPointShootingLowThrustPhase(const std::string& name,
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

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount) = 0;

            virtual void output_STMs();

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
            //calcbounds
            //these (can) go in the specialized phase
            virtual void calcbounds_phase_main() = 0;

            //these happen only at the TwoPointShootingLowThrustPhase level
            void calcbounds_control();
            void calcbounds_virtual_propellant_tanks();

            void calcbounds_deltav_contribution();

                        
            virtual void calcbounds_match_point_constraints(); //TwoPointShootingLowThrustPhase can do the sparsity pattern except for Sundman phases
            virtual void calcbounds_distance_constraints(); //TwoPointShootingLowThrustPhase can do the sparsity pattern except for Sundman phases

            //process
            void process_control(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //these go in the specialized phase
            virtual void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;//TwoPointShootingLowThrustPhase can't handle these derivatives on its own - specialized class has to do it

            virtual void process_distance_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;//TwoPointShootingLowThrustPhase can't handle these derivatives on its own - specialized class has to do it

            //fields
            size_t num_timesteps;
            std::vector<doubleType> ThrustStepLengths;
            std::vector<double> dThrustStepLengths_dPropagationVariable;
            std::vector<doubleType> ForwardPropagationStepTimes;
            std::vector<doubleType> BackwardPropagationStepTimes;
            std::vector<double> dForwardPropagationStepTimes_dPropagationVariable;
            std::vector<double> dBackwardPropagationStepTimes_dPropagationVariable;

            //control
            std::vector< math::Matrix<doubleType> > ControlVector;
            std::vector<doubleType> throttle;
            size_t num_controls; //3 normally, but 4 for VSI systems
            bool isVSI;
            std::vector<size_t> TruthTable_MatchConstraints_Derivative_wrt_Control;

            //thruster performance
            std::vector<doubleType> max_thrust;
            std::vector<doubleType> max_mass_flow_rate;
            std::vector<doubleType> Isp;
            std::vector<doubleType> power;
            std::vector<doubleType> active_power;
            std::vector<size_t> number_of_active_engines;
            std::vector<int> ThrottleLevel;
            std::vector<std::string> ThrottleLevelString;

            //delta-v
            std::vector<doubleType> stepDeltavMagnitude;
            std::vector< math::Matrix<doubleType> > stepDeltaV;

            //forced coasts
            math::Matrix<doubleType> state_after_initial_coast;
            math::Matrix<double> STM_initial_coast;
            math::Matrix<double> STM_Augmented_initial_coast;
            math::Matrix<double> InitialCoast_dStatedIndependentVariable;
            math::Matrix<doubleType> state_before_terminal_coast;
            math::Matrix<double> STM_terminal_coast;
            math::Matrix<double> STM_Augmented_terminal_coast;
            math::Matrix<double> TerminalCoast_dStatedIndependentVariable;

            //distance constraint
            size_t number_of_distance_constraints;
            std::vector < std::tuple<int, double, double> > distance_constraint_definitions; //body, lower bound in km, upper bound in km
            std::vector < std::vector < math::Matrix<doubleType> > > distance_constraint_relative_position; //step, constraint, state variable
            std::vector < std::vector < math::Matrix<double> > > distance_constraint_body_position_time_derivatives;
            std::vector< std::vector<doubleType> > distance_from_body; //step, constraint

            //derivative indices
            std::vector< std::vector<size_t> > Xindex_control;//stepIndex, controlIndex
            std::vector< std::vector<size_t> > Gindices_fixed_inertial_control_wrt_this_step_control;//step, control variable (num_steps - 1 long)
            std::vector< std::vector<size_t> > Gindices_fixed_inertial_control_wrt_first_step_control;//step, control variable (num_steps - 1 long)
            std::vector< std::vector< std::vector<size_t> > > G_indices_match_point_constraints_wrt_Control;//state, step, control variable

            std::vector < std::vector<size_t> > G_indices_control_magnitude_constraints_wrt_Control;

            std::vector< std::vector< std::vector< std::vector<size_t> > > > G_indices_distance_constraints_wrt_ForwardControl; //step, constraint, controlstep, control variable
            std::vector< std::vector< std::vector< std::vector<size_t> > > > G_indices_distance_constraints_wrt_BackwardControl; //step, constraint, controlstep, control variable
            std::vector< std::vector< std::vector<size_t> > > G_indices_distance_constraints_wrt_LeftBoundaryState; //step, constraint, variable
            std::vector< std::vector< std::vector<size_t> > > G_indices_distance_constraints_wrt_RightBoundaryState; //step, constraint, variable
            std::vector< std::vector< std::vector<size_t> > > G_indices_distance_constraints_wrt_LeftBoundaryTime; //step, constraint, variable
            std::vector< std::vector< std::vector<size_t> > > G_indices_distance_constraints_wrt_RightBoundaryTime; //step, constraint, variable
            std::vector< std::vector< size_t> > G_indices_distance_constraints_wrt_PhaseFlightTime;

            std::vector<size_t> dIndex_virtual_electric_propellant_StateAfterDepartureDecisionVariables;
            std::vector<size_t> dIndex_virtual_electric_propellant_StateBeforeArrivalDecisionVariables;
            std::vector<size_t> dIndex_virtual_electric_propellant_StateAfterDepartureDecisionVariables_wrt_Time;
            std::vector<size_t> dIndex_virtual_electric_propellant_StateBeforeArrivalDecisionVariables_wrt_Time;
            std::vector<size_t> Gindex_virtual_electric_propellant_StateAfterDepartureDecisionVariables;
            std::vector<size_t> Gindex_virtual_electric_propellant_StateBeforeArrivalDecisionVariables;
            std::vector<size_t> Gindex_virtual_electric_propellant_StateAfterDepartureDecisionVariables_wrt_Time;
            std::vector<size_t> Gindex_virtual_electric_propellant_StateBeforeArrivalDecisionVariables_wrt_Time;

            std::vector<size_t> G_indices_deltav_wrt_LeftBoundaryState;
            std::vector<size_t> G_indices_deltav_wrt_RightBoundaryState;
            std::vector<size_t> G_indices_deltav_wrt_LeftBoundaryState_wrt_Time;
            std::vector<size_t> G_indices_deltav_wrt_RightBoundaryState_wrt_Time;
            std::vector< std::vector<size_t> > G_indices_deltav_wrt_Control;
            size_t G_indices_deltav_wrt_PhaseFlightTime;
            
        };
    }//close namespace Phases
}//close namespace EMTG