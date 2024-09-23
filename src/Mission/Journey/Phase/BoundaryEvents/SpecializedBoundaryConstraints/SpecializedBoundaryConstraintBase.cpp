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

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            SpecializedBoundaryConstraintBase::SpecializedBoundaryConstraintBase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition ) :
                SpecializedBoundaryConstraintBase()
            {

                this->name = name;
                this->journeyIndex = journeyIndex;
                this->phaseIndex = phaseIndex;
                this->stageIndex = stageIndex;
                this->myUniverse = Universe;
                this->mySpacecraft = mySpacecraft;
                this->myOptions = myOptions;
                this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];
                this->myBoundaryEvent = myBoundaryEvent;
                this->constraintDefinition = constraintDefinition;
                this->myReferenceFrame = ReferenceFrame::ICRF;
                this->ReferenceFrameStrings = std::vector<std::string> ({ "ICRF", "J2000_BCI", "J2000_BCF", "TrueOfDateBCI", "TrueOfDate_BCF", "PrincipleAxes", "Topocentric", "Polar" });
            }

            void SpecializedBoundaryConstraintBase::setup_calcbounds(
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
            }
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG