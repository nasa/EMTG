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

//MGAnDSMs maneuver epoch constraint
//9-22-2017

#pragma once

#include "doubleType.h"
#include "sparsey_thing.h"

#include "universe.h"
#include "Spacecraft.h"

//#include "MGAnDSMs_subphase.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare parent class
        class MGAnDSMs_subphase;

        class MGAnDSMs_maneuver_constraint : public sparsey_thing
        {
        public:
            //constructor
            MGAnDSMs_maneuver_constraint() {};
            MGAnDSMs_maneuver_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_subphase* mySubPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

			virtual ~MGAnDSMs_maneuver_constraint();

            //clone
            virtual MGAnDSMs_maneuver_constraint* clone() const = 0;

            //calcbounds goes in the specialized phase
            virtual void calcbounds() = 0;

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

            //process goes in the specialized phase
            virtual void process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

        protected:

            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            size_t subphaseIndex;
            MGAnDSMs_subphase* mySubPhase;
            Astrodynamics::universe* Universe;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;
            JourneyOptions* myJourneyOptions;
            std::string ConstraintDefinition;

            double lowerBound, upperBound;
        };
        
        inline MGAnDSMs_maneuver_constraint * new_clone(MGAnDSMs_maneuver_constraint const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG