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

//EMTGv9 TwoPointShootingPhase
//Jacob Englander 6-22-2017

#pragma once

#include "phase.h"

namespace EMTG 
{
    namespace Phases
    {
        class TwoPointShootingPhase : public phase
        {
        public:
            //constructor
            TwoPointShootingPhase() : phase::phase() {};

            TwoPointShootingPhase(const std::string& name,
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

            //these go in the specialized phase
            virtual void calcbounds_phase_main() = 0;
            virtual void calcbounds_match_point_constraints() = 0;

            void process_phase_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_forward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_backward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;//specialized class forms STMs, MTMs, and HPTMs, base class does the boundary derivative only

            //fields
            math::Matrix<doubleType> match_point_state_minus;
            math::Matrix<doubleType> match_point_state_plus;
            std::vector< math::Matrix<doubleType> > spacecraft_state_event_minus;
            std::vector< math::Matrix<doubleType> > spacecraft_state_event_plus;

            math::Matrix<double> ForwardHPTM;
            math::Matrix<double> BackwardHPTM;

            std::vector< math::Matrix<double> > ForwardSTM;
            std::vector< math::Matrix<double> > BackwardSTM;
            std::vector< math::Matrix<double> > Forward_dStatedIndependentVariable;
            std::vector< math::Matrix<double> > Backward_dStatedIndependentVariable;
            std::vector< math::Matrix<double> > ForwardAugmentedSTM; //includes mass, time, control if necessary, etc
            std::vector< math::Matrix<double> > BackwardAugmentedSTM; //includes mass, time, control if necessary, etc
        
            //constraint indices
            std::vector<size_t> Findices_match_point_constraints;

            //derivative indices
            std::vector<size_t> ListOfVariablesAffectingLeftBoundary;
            std::vector<size_t> ListOfVariablesAffectingRightBoundary;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryVariables;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryVariables;
            std::vector<size_t> ListOfTimeVariablesAffectingLeftBoundary;
            std::vector<size_t> ListOfTimeVariablesAffectingRightBoundary;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByTimeVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByTimeVariable; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables;
            std::vector<size_t> G_indices_match_point_constraints_wrt_PhaseFlightTime;

            std::vector< std::vector<bool> > LeftBoundaryVariableAffectsMatchConstraint;//constraintIndex, entry in ListOfVariablesAffectingLeftBoundary
            std::vector< std::vector<bool> > RightBoundaryVariableAffectsMatchConstraint;//constraintIndex, entry in ListOfVariablesAffectingRightBoundary
            std::vector< std::vector<bool> > LeftBoundaryTimeVariableAffectsMatchConstraint;//constraintIndex, entry in ListOfTimeVariablesAffectingLeftBoundary
            std::vector< std::vector<bool> > RightBoundaryTimeVariableAffectsMatchConstraint;//constraintIndex, entry in ListOfTimeVariablesAffectingRightBoundary
            
        };
    }//close namespace Phases
}//close namespace EMTG