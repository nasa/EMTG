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

//parallel shooting finite burn (PSFB) step for EMTGv9
//Jacob Englander 2-23-2018

#pragma once

#include "doubleType.h"

#include "ParallelShootingStep.h"

#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare PSFBphase
        //class PSFBphase;

        class PSFBstep : virtual public ParallelShootingStep
        {
        public:
            //constructor
            PSFBstep();
            PSFBstep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual PSFBstep* clone() const { return new PSFBstep(*this); }

            //destructor
            virtual ~PSFBstep();

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            virtual void calcbounds_step();

            virtual void process_step(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //configure propagator
            virtual void configure_propagator();

            //calcbounds
            virtual void calcbounds_deltav_contribution() {}; //empty for now

            //process

            virtual void process_step_main(const std::vector<doubleType>& X,
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
                const bool& needG) {};//empty for now

            //fields
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;  
            size_t total_number_of_states_to_integrate;
        };

        inline PSFBstep * new_clone(PSFBstep const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG