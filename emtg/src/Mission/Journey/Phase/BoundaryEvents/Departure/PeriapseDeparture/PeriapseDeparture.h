
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

#include "arrival.h"
#include "departure.h"
#include "PeriapseBoundary.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class PeriapseDeparture : virtual public DepartureEvent,
            virtual public PeriapseBoundary
        {
        public:
            //default constructor
            PeriapseDeparture() {};

            //specialized constructor
            PeriapseDeparture(const std::string& name,
                           const size_t& journeyIndex,
                           const size_t& phaseIndex,
                           size_t& stageIndex,
                           Astrodynamics::universe* Universe,
                           HardwareModels::Spacecraft* mySpacecraft,
                           missionoptions* myOptions,
                           ArrivalEvent* PreviousPhaseArrivalEvent);

            virtual void initialize(const std::string& name,
                                    const size_t& journeyIndex,
                                    const size_t& phaseIndex,
                                    size_t& stageIndex,
                                    Astrodynamics::universe* Universe,
                                    HardwareModels::Spacecraft* mySpacecraft,
                                    missionoptions* myOptions,
                                    ArrivalEvent* PreviousPhaseArrivalEvent);

            //destructor
            virtual ~PeriapseDeparture() {};

            //output
            virtual void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount) = 0;

        protected:

            //method to calculate event left side

            virtual void calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                    const std::vector<double>& VelocityMagnitudeBounds,
                                                    std::vector<size_t> timeVariables);
            
            virtual void calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                    const std::vector<double>& VelocityMagnitudeBounds,
                                                    const std::vector<double>& MassBounds,
                                                    std::vector<size_t> timeVariables);

            virtual void calcbounds_event_right_side();

            virtual void calcbounds_specialized_constraints();

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_event_right_side(const std::vector<doubleType>& X, //this is just to handle the derivative with respect to wait time
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG