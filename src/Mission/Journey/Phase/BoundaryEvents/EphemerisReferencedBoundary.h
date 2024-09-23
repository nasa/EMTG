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

#include "BoundaryEventBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisReferencedBoundary : public virtual BoundaryEventBase
        {
        public:
            //constructor
            EphemerisReferencedBoundary() : BoundaryEventBase::BoundaryEventBase() {};
            EphemerisReferencedBoundary(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                                    const size_t& journeyIndex,
                                    const size_t& phaseIndex,
                                    size_t& stageIndex,
                                    Astrodynamics::universe* Universe,
                                    HardwareModels::Spacecraft* mySpacecraft,
                                    missionoptions* myOptions);

        protected:
            //these go in the specialized event
            virtual void calcbounds_event_interface_state(const std::vector<double>& RAbounds,
                const std::vector<double>& DECbounds,
                std::vector<double>& MassBounds,
                const std::vector<double>& EpochBounds,
                std::vector<size_t>& timeVariables); //I did not make this abstract because it gets overloaded by two of the three derived classes

            virtual void process_event_interface_state(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //fields
            double semi_axis_a;
            double semi_axis_b;
            double semi_axis_c;

            math::Matrix<doubleType> state_on_interface_spherical;
            math::Matrix<doubleType> state_on_interface_cartesian;

            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_on_interface_cartesian;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_state_on_interface_cartesian_wrt_Time;//Xindex, stateIndex, derivative value

            std::vector<size_t> dIndex_interface_position_wrt_RA;//x, y, z
            std::vector<size_t> dIndex_interface_position_wrt_DEC;//x, y, z
            
        };
    }//close namespace events
}//close namespace EMTG