
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
        class PeriapseBoundary : virtual public BoundaryEventBase
        {
        public:
            //default constructor
            PeriapseBoundary();

            //specialized constructor
            PeriapseBoundary(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //initialize method
            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

        protected:
            //these go in the specialized event
            virtual void calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                const std::vector<double>& VelocityMagnitudeBounds,
                const std::vector<double>& MassBounds);

            virtual void calcbounds_event_right_side() = 0;

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //fields

            size_t dIndex_mass_wrt_mass;
            size_t Xindex_rMag;
            size_t Xindex_RA;
            size_t Xindex_DEC;
            size_t Xindex_vMag;
            size_t Xindex_vRA;
            size_t Xindex_vDEC;
            size_t Xindex_AZ;
            size_t Xindex_FPA;
            size_t Xindex_mass;

            std::vector<size_t> dIndex_periapse_state_wrt_r;//x, y, z
            std::vector<size_t> dIndex_periapse_state_wrt_RA;//x, y, xdot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_DEC;//x, y, z, xdot, ydot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_v;//xdot, ydot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_vRA;//xdot, ydot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_vDEC;//xdot, ydot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_AZ;//xdot, ydot, zdot
            std::vector<size_t> dIndex_periapse_state_wrt_FPA;//xdot, ydot, zdot

            size_t Gindex_rdotv_wrt_RA;
            size_t Gindex_rdotv_wrt_DEC;
            size_t Gindex_rdotv_wrt_vRA;
            size_t Gindex_rdotv_wrt_vDEC;
            size_t Gindex_rdotv_wrt_AZ;
            size_t Gindex_rdotv_wrt_FPA;
        };//end class PeriapseBoundary

    }//end namespace BoundaryEvents
}//end namespace EMTG