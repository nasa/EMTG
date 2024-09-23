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

//EMTGv9 CoastPhase
//Jacob Englander 11-20-2017

#pragma once

#include "TwoPointShootingPhase.h"
#include "IntegratedPropagator.h"
#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG
{
    namespace Phases
    {
        class CoastPhase : public TwoPointShootingPhase
        {
        public:
            //constructor
            CoastPhase() : TwoPointShootingPhase::TwoPointShootingPhase() {};
            CoastPhase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            CoastPhase(const std::string& name,
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

            virtual ~CoastPhase();

            virtual void initialize();

            //clone
            virtual CoastPhase* clone() const { return new CoastPhase(*this); }

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_STMs();

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            //calcbounds
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

            void calcbounds();

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //calcbounds
            void calcbounds_phase_main();

            virtual void calcbounds_match_point_constraints();

            virtual void calcbounds_virtual_propellant_tanks(); //ACS

            void calcbounds_deltav_contribution();

            //process
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
                const bool& needG);

            virtual void process_backward_half_phase(const std::vector<doubleType>& X,
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
                const bool& needG);

            void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);//ACS

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields - aka propagators
            size_t num_timesteps;
            

            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            
            Integration::IntegrationScheme* myIntegrationScheme;

            Astrodynamics::PropagatorBase* ForwardHalfPhasePropagator;
            Astrodynamics::PropagatorBase* BackwardHalfPhasePropagator;
            double ForwardIntegrationStepLength, BackwardIntegrationStepLength;
            double MatchPointFraction;

            size_t total_number_of_states_to_integrate;

            math::Matrix<double> ForwardSTM, BackwardSTM, ForwardSPTM, BackwardSPTM;
            math::Matrix<double> Forward_dPropagatedStatedIndependentVariable, Backward_dPropagatedStatedIndependentVariable;
            double dForwardStepSize_dPropagationVariable;
            double dBackwardStepSize_dPropagationVariable;

			Astrodynamics::TimeDomainSpacecraftEOM myEOM;

            bool isKeplerian;
        };
    }//close namespace Phases
}//close namespace EMTG