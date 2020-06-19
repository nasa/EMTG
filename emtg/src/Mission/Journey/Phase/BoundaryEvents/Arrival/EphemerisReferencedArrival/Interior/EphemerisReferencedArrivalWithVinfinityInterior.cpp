
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

#include "EphemerisReferencedArrivalWithVinfinityInterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedArrivalWithVinfinityInterior::EphemerisReferencedArrivalWithVinfinityInterior(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->EphemerisReferencedArrivalInterior::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            if (this->isLastEventInJourney)
            {
                if (this->myOptions->Journeys[this->journeyIndex].journey_end_TCM > math::SMALL)
                {
                    this->hasTCM = true;
                    this->TCM_magnitude = this->myOptions->Journeys[this->journeyIndex].journey_end_TCM;
                }
            }
            else if (this->myOptions->TCM_pre_flyby > math::SMALL)
            {
                this->hasTCM = true;
                this->TCM_magnitude = this->myOptions->TCM_pre_flyby;
            }
            else
                this->hasTCM = false;
        }

        //******************************************calcbounds methods

        void EphemerisReferencedArrivalWithVinfinityInterior::calcbounds_event_interface_state(const std::vector<double>& vinf_max,
            const std::vector<double>& RAbounds,
            const std::vector<double>& DECbounds,
            std::vector<double>& MassBounds,
            const std::vector<double>& EpochBounds,
            std::vector<size_t> timeVariables)
        {
            //Step 1: calcbounds for the v-infinity, also in polar coordinates
            this->Xlowerbounds->push_back(vinf_max[0]);
            this->Xupperbounds->push_back(vinf_max[1]);
            this->Xdescriptions->push_back(this->prefix + "event interface state vMAG");
            for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                this->Derivatives_of_state_on_interface_cartesian.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                this->dIndex_interface_velocity_wrt_v.push_back(this->Derivatives_of_state_on_interface_cartesian.size() - 1);
            }
            X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);

            this->Xlowerbounds->push_back(-8.0 * math::PI);
            this->Xupperbounds->push_back(8.0 * math::PI);
            this->Xdescriptions->push_back(this->prefix + "event interface state vRA");
            for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                this->Derivatives_of_state_on_interface_cartesian.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                this->dIndex_interface_velocity_wrt_vRA.push_back(this->Derivatives_of_state_on_interface_cartesian.size() - 1);
            }
            this->X_scale_factors->push_back(1.0);

            this->Xlowerbounds->push_back(DECbounds[0] * math::PI);
            this->Xupperbounds->push_back(DECbounds[1] * math::PI);
            this->Xdescriptions->push_back(this->prefix + "event interface state vDEC");
            for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                this->Derivatives_of_state_on_interface_cartesian.push_back({ this->Xdescriptions->size() - 1, stateIndex, 1.0 });
                this->dIndex_interface_velocity_wrt_vDEC.push_back(this->Derivatives_of_state_on_interface_cartesian.size() - 1);
            }
            this->X_scale_factors->push_back(1.0);

            //Step 2: call the base class
            this->EphemerisReferencedArrivalInterior::calcbounds_event_interface_state(RAbounds, DECbounds, MassBounds, EpochBounds, timeVariables);
        }//end calcbounds_event_interface_state()


        void EphemerisReferencedArrivalWithVinfinityInterior::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //base class
            this->EphemerisReferencedArrivalInterior::calcbounds_event_left_side(timeVariables);

            //no specialized variables
        }//end calcbounds_event_left_side()

        void EphemerisReferencedArrivalWithVinfinityInterior::calcbounds_event_right_side()
        {
            //base class
            this->EphemerisReferencedArrivalInterior::calcbounds_event_right_side();

            //no specialized variables
        }//end calcbounds_event_right_side()

         //******************************************process methods

        void EphemerisReferencedArrivalWithVinfinityInterior::process_event_interface_state(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: velocity
            this->state_on_interface_spherical(3) = X[Xindex++];
            this->state_on_interface_spherical(4) = X[Xindex++];
            this->state_on_interface_spherical(5) = X[Xindex++];

            //Step 2: convert velocity to cartesian
            doubleType& v = this->state_on_interface_spherical(3);
            doubleType& vRA = this->state_on_interface_spherical(4);
            doubleType& vDEC = this->state_on_interface_spherical(5);
            doubleType cosvRA = cos(vRA);
            doubleType sinvRA = sin(vRA);
            doubleType cosvDEC = cos(vDEC);
            doubleType sinvDEC = sin(vDEC);

            this->state_on_interface_cartesian(3) = v * cosvRA * cosvDEC;
            this->state_on_interface_cartesian(4) = v * sinvRA * cosvDEC;
            this->state_on_interface_cartesian(5) = v * sinvDEC;

            //Step 3: derivatives of velocity
            if (needG)
            {
                double dvx_dv = (cosvRA * cosvDEC) _GETVALUE;
                double dvx_dvRA = (-v * sinvRA * cosvDEC)_GETVALUE;
                double dvx_dvDEC = (-v * cosvRA * sinvDEC)_GETVALUE;
                double dvy_dv = (sinvRA * cosvDEC) _GETVALUE;
                double dvy_dvRA = (v * cosvRA * cosvDEC)_GETVALUE;
                double dvy_dvDEC = (-v * sinvRA * sinvDEC)_GETVALUE;
                double dvz_dv = sinvDEC _GETVALUE;
                double dvz_dvRA = 0.0;
                double dvz_dvDEC = (v * cosvDEC) _GETVALUE;

                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_v[0]]) = dvx_dv;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_v[1]]) = dvy_dv;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_v[2]]) = dvz_dv;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vRA[0]]) = dvx_dvRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vRA[1]]) = dvy_dvRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vRA[2]]) = dvz_dvRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vDEC[0]]) = dvx_dvDEC;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vDEC[1]]) = dvy_dvDEC;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_velocity_wrt_vDEC[2]]) = dvz_dvDEC;
            }

            //Step 4: position (base class)
            this->EphemerisReferencedArrivalInterior::process_event_interface_state(X, Xindex, F, Findex, G, needG);
        }//end process_event_interface_state()

    }//end namespace BoundaryEvents
}//end namespace EMTG