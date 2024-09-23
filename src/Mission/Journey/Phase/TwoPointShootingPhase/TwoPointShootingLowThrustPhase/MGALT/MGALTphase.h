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

//EMTGv9 MGALTphase
//Jacob Englander 6-24-2017

#pragma once

#include "TwoPointShootingLowThrustPhase.h"

#include "../Propagation/KeplerPropagatorTimeDomain.h"
#include "ForwardBoundedImpulseManeuver.h"
#include "BackwardBoundedImpulseManeuver.h"

namespace EMTG
{
    namespace Phases
    {
        class MGALTphase : public TwoPointShootingLowThrustPhase
        {
        public:
            //constructor
            MGALTphase() : TwoPointShootingLowThrustPhase::TwoPointShootingLowThrustPhase() {};
            MGALTphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            //destructor
            virtual ~MGALTphase();

            //clone
            virtual MGALTphase* clone() const { return new MGALTphase(*this); }

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
            Astrodynamics::KeplerPropagatorTimeDomain InitialCoastPropagatorObject;
            Astrodynamics::KeplerPropagatorTimeDomain TerminalCoastPropagatorObject;
            std::vector< Astrodynamics::KeplerPropagatorTimeDomain > ForwardPropagatorObjects;
            std::vector< Astrodynamics::KeplerPropagatorTimeDomain > BackwardPropagatorObjects;

            std::vector< ForwardBoundedImpulseManeuver > ForwardManeuverObjects;
            std::vector< BackwardBoundedImpulseManeuver > BackwardManeuverObjects;

            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;

            std::vector< math::Matrix<double> > ForwardMTM;
            std::vector< math::Matrix<double> > BackwardMTM;
            std::vector < math::Matrix<double> > MTM_STM_blocks_Forward;
            std::vector < math::Matrix<double> > MTM_STM_blocks_Backward;
            std::vector < math::Matrix<double> > MTM_STM_cumulative_chain_Forward;
            std::vector < math::Matrix<double> > MTM_STM_cumulative_chain_Backward;
            
            std::vector<double> dImpulseEpoch_dPropagationVariable;

            //distance constraint
            std::vector < std::vector<math::Matrix<double> > > Forward_cumulative_STM_MTM_triangle;
            std::vector < std::vector<math::Matrix<double> > > Backward_cumulative_STM_MTM_triangle;
            std::vector< math::Matrix<double> > Forward_boundary_STM_MTM;
            std::vector< math::Matrix<double> > Backward_boundary_STM_MTM;

        };
    }//close namespace Phases
}//close namespace EMTG