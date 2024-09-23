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

#include "ObjectiveFunctionBase.h"

namespace EMTG
{
    namespace ObjectiveFunctions
    {
        class ArriveAsEarlyAsPossibleObjective : public ObjectiveFunctionBase
        {
        public:
            //constructors
            ArriveAsEarlyAsPossibleObjective() {};

            ArriveAsEarlyAsPossibleObjective(Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                Mission* myMission);

            //destructor
            virtual ~ArriveAsEarlyAsPossibleObjective() {};


            //public methods
            void calcbounds();

            void process(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void output(std::ofstream& outputfile);

        protected:
            //fields
            std::vector<size_t> Xindices_ArrivalEpoch;

            BoundaryEvents::ArrivalEvent* myArrivalEvent;                
        };
    }//end namespace ObjectiveFunctions
}//end namespace EMTG