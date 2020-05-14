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


//ParallelShootingStep maneuver constraint base class
//Jacob Englander 4-13-2018

#include "ParallelShootingStep_maneuver_constraint.h"

#include "boost/algorithm/string/split.hpp"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStep_maneuver_constraint::ParallelShootingStep_maneuver_constraint(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& subStepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* myStep,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            ParallelShootingStep_maneuver_constraint()
        {
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stepIndex = stepIndex;
            this->subStepIndex = subStepIndex;
            this->stageIndex = stageIndex;
            this->myUniverse = myUniverse;
            this->mySpacecraft = mySpacecraft;
            this->myOptions = myOptions;
            this->myJourneyOptions = &(this->myOptions->Journeys[this->journeyIndex]);
            this->myStep = myStep;
            this->ConstraintDefinition = ConstraintDefinition;
            this->name = this->myStep->getName();
        }//end constructor

		ParallelShootingStep_maneuver_constraint::~ParallelShootingStep_maneuver_constraint() {}

        void ParallelShootingStep_maneuver_constraint::setup_calcbounds(std::vector<double>* Xupperbounds,
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
            std::vector<double>* A)
        {
            this->prefix = this->name + ": ";

            this->sparsey_thing::setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
				F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);
        }//end setup_calcbounds()
    }//close namespace Phases
}//close namespace EMTG