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

//EMTGv9 FBLTphase
//Jacob Englander 1-2-2018

#pragma once

#include "TwoPointShootingLowThrustPhase.h"

#include "IntegratedPropagator.h"
#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG
{
    namespace Phases
    {
        class FBLTphase : public TwoPointShootingLowThrustPhase
        {
        public:
            //constructor
            FBLTphase() : TwoPointShootingLowThrustPhase::TwoPointShootingLowThrustPhase() {};
            FBLTphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            virtual ~FBLTphase();

            //clone
            virtual FBLTphase* clone() const { return new FBLTphase(*this); }

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            //calcbounds goes in the specialized phase
            virtual void calcbounds();

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //calcbounds

            //these (can) go in the specialized phase
            virtual void calcbounds_phase_main();

            //process

            //these go in the specialized constraint but they also exist for MGALT
            virtual void process_phase_flight_time(const std::vector<doubleType>& X,
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

            virtual void process_distance_constraints(const std::vector<doubleType>& X,
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

            //fields
            Astrodynamics::PropagatorBase* InitialCoastPropagatorObject;
            Astrodynamics::PropagatorBase* TerminalCoastPropagatorObject;
            boost::ptr_vector < Astrodynamics::PropagatorBase > ForwardPropagatorObjects;
            boost::ptr_vector < Astrodynamics::PropagatorBase > BackwardPropagatorObjects;
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;
            double ForcedCoast_dStepSize_dPropagationVariable;

            size_t total_number_of_states_to_integrate;

            std::vector < math::Matrix<double> > ForwardStrippedAugmentedSTM;
            std::vector < math::Matrix<double> > BackwardStrippedAugmentedSTM;
            std::vector < math::Matrix<double> > STM_cumulative_chain_Forward;
            std::vector < math::Matrix<double> > STM_cumulative_chain_Backward;
            std::vector < math::Matrix<double> > STM_stripped_cumulative_chain_Forward;
            std::vector < math::Matrix<double> > STM_stripped_cumulative_chain_Backward;

            //distance constraint
            std::vector < std::vector<math::Matrix<double> > > Forward_cumulative_STM_triangle;
            std::vector < std::vector<math::Matrix<double> > > Backward_cumulative_STM_triangle;
            std::vector < std::vector<math::Matrix<double> > > Forward_cumulative_stripped_STM_triangle;
            std::vector < std::vector<math::Matrix<double> > > Backward_cumulative_stripped_STM_triangle;
            std::vector < math::Matrix<double> > Forward_boundary_STM;
            std::vector < math::Matrix<double> > Backward_boundary_STM;

        };
    }//close namespace Phases
}//close namespace EMTG