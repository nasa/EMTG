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

#pragma once

#include "doubleType.h"

#include <string>

#include "universe.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "sparsey_thing.h"

#include "BoundaryEventBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //forward declare stuff
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {

            class SpecializedBoundaryConstraintBase : public sparsey_thing
            {
            public:
                //constructors
                SpecializedBoundaryConstraintBase() {};

                SpecializedBoundaryConstraintBase(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~SpecializedBoundaryConstraintBase() {};


                //public methods

                virtual void output(std::ofstream& outputfile) = 0;

                //calcbounds
                void setup_calcbounds(
                    std::vector<double>* Xupperbounds,
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

                //calcbounds goes in the specialized event
                virtual void calcbounds() = 0;

                //process
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
                Astrodynamics::universe* myUniverse;
                Astrodynamics::body* myBody;
                HardwareModels::Spacecraft* mySpacecraft;
                missionoptions* myOptions;
                JourneyOptions* myJourneyOptions;
                BoundaryEventBase* myBoundaryEvent;
                std::string constraintDefinition;    

                ReferenceFrame myReferenceFrame;
                std::vector<std::string> ReferenceFrameStrings;
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG