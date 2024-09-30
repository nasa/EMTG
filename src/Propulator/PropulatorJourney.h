// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

//PropulatorJourney

#include "missionoptions.h"
#include "journeyoptions.h"
#include "universe.h"
#include "Spacecraft.h"

#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"
#include "PropagatorBase.h"

namespace EMTG
{
    namespace Propulator
    {
        class PropulatorJourney
        {
        public:
            PropulatorJourney();
            PropulatorJourney(const size_t& journeyIndex,
                size_t& stageIndex,
                missionoptions* options,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* spacecraft);

            ~PropulatorJourney();

            void initialize(const size_t& journeyIndex,
                size_t& stageIndex,
                missionoptions* options,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* spacecraft);

            void setStepSize(const double& stepSize) { this->myPropagator->setPropagationStepSize(stepSize); }

            math::Matrix<double> getSTM() { return this->STM; }

            void propulate(const math::Matrix<double>& InputState,
                const math::Matrix<double>& ControlVector,
                const double& DutyCycle,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

            void propulate(const math::Matrix<double>& InputState,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

        protected:
            //structure
            size_t journeyIndex;
            size_t stageIndex;
            JourneyOptions* myJourneyOptions;
            Astrodynamics::universe* myUniverse;
            missionoptions* myOptions;
            HardwareModels::Spacecraft* mySpacecraft;

            //states
            math::Matrix<double> StateBeforePropagation;
            math::Matrix<double> StateAfterPropagation;
            math::Matrix<double> STM;
            math::Matrix<double> dPropagatedStatedIndependentVariable;

            //propagator and dynamics
            Astrodynamics::PropagatorBase* myPropagator;
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;
            size_t total_number_of_states_to_integrate;
            double dPropagationTime_dPropagationTime;
        };
    }//close namespace Propulator
}//close namespace EMTG