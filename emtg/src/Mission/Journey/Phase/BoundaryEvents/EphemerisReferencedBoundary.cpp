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

#include "EphemerisReferencedBoundary.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedBoundary::EphemerisReferencedBoundary(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            state_on_interface_spherical(8, 1, 0.0),
            state_on_interface_cartesian(8, 1, 0.0)
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions);
        }
        
        void EphemerisReferencedBoundary::initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions)
        {
            //base class
            this->BoundaryEventBase::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            state_on_interface_spherical.resize(8, 1, 0.0);
            state_on_interface_cartesian.resize(8, 1, 0.0);
        }//end initialize()

        //calcbounds methods

        void EphemerisReferencedBoundary::calcbounds_event_interface_state(const std::vector<double>& RAbounds,
            const std::vector<double>& DECbounds, 
            std::vector<double>& MassBounds, 
            const std::vector<double>& EpochBounds,
            std::vector<size_t>& timeVariables)
        {
            //Step 1: What is the first decision variable in this event?
            for (int Xindex = this->Xdescriptions->size() - 1; Xindex >= 0; --Xindex)
            {
                if (this->Xdescriptions->at(Xindex).find(this->prefix) < 1024)
                {
                    this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();
                    break;
                }
            }

            //Step 2: encode position (assume that the calling function has already encoded spherical velocity)
            //this is just two angles in polar coordinates - RA measured in the x-y plane from x and DEC measured up from the x-y plane to the +z axis
            //yes, there is a singularity at the pole but we can live with that
            this->Xlowerbounds->push_back(RAbounds[0] * math::deg2rad);
            this->Xupperbounds->push_back(RAbounds[1] * math::deg2rad);
            this->Xdescriptions->push_back(this->prefix + "event interface state RA");
            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                this->Derivatives_of_state_on_interface_cartesian.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                this->dIndex_interface_position_wrt_RA.push_back(this->Derivatives_of_state_on_interface_cartesian.size() - 1);
            }
            this->X_scale_factors->push_back(1.0);

            this->Xlowerbounds->push_back(DECbounds[0] * math::deg2rad);
            this->Xupperbounds->push_back(DECbounds[1] * math::deg2rad);
            this->Xdescriptions->push_back(this->prefix + "event interface state DEC");
            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                this->Derivatives_of_state_on_interface_cartesian.push_back({ this->Xdescriptions->size() - 1, stateIndex, 1.0 });
                this->dIndex_interface_position_wrt_DEC.push_back(this->Derivatives_of_state_on_interface_cartesian.size() - 1);
            }
            this->X_scale_factors->push_back(1.0);

            //Step 3: mass variable
            this->Xlowerbounds->push_back(MassBounds[0]);
            this->Xupperbounds->push_back(MassBounds[1]);
            this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
            this->Xdescriptions->push_back(this->prefix + "event left state mass");
            size_t Xindex_mass = this->Xdescriptions->size() - 1;
            this->Derivatives_of_state_on_interface_cartesian.push_back({ Xindex_mass, 6, 1.0 });
            this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_state_on_interface_cartesian.size() - 1;


            //Step 4: left epoch
            if (this->isFirstEventInMission)
            {
                this->Xlowerbounds->push_back(EpochBounds[0]);
                this->Xupperbounds->push_back(EpochBounds[1]);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(this->prefix + "event left state epoch");
                this->Xindices_EventLeftEpoch.push_back(this->Xdescriptions->size() - 1);
                this->Derivatives_of_state_on_interface_cartesian_wrt_Time.push_back({ this->Xdescriptions->size() - 1, 7, 1.0 });
            }

            //the epoch of the interface state depends on all previous epochs - this is stored slightly differently compared to other boundary types
            for (size_t Xindex : timeVariables)
            {
                this->Xindices_EventLeftEpoch.push_back(Xindex);
                this->Derivatives_of_state_on_interface_cartesian_wrt_Time.push_back({ Xindex, 7, 1.0 });
            }
        }//end calcbounds_event_left_side

        void EphemerisReferencedBoundary::process_event_interface_state(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract RA and DEC, form the position vector
            //we will assume that whatever function is calling this one has already set the ephemeris-referenced velocity components!
            doubleType RA = X[Xindex++];
            doubleType DEC = X[Xindex++];

            doubleType cosRA = cos(RA);
            doubleType sinRA = sin(RA);
            doubleType cosDEC = cos(DEC);
            doubleType sinDEC = sin(DEC);
            doubleType cosRA2 = cosRA * cosRA;
            doubleType sinRA2 = sinRA * sinRA;
            doubleType cosDEC2 = cosDEC * cosDEC;
            doubleType sinDEC2 = sinDEC * sinDEC;
            double a2 = this->semi_axis_a * this->semi_axis_a;
            double b2 = this->semi_axis_b * this->semi_axis_b;
            double c2 = this->semi_axis_c * this->semi_axis_c;

            doubleType r = sqrt(1.0 / ( (cosRA2 * cosDEC2 / a2)
                                      + (sinRA2 * cosDEC2 / b2)
                                      + (sinDEC2 / c2) ) );

            this->state_on_interface_spherical(0) = r;
            this->state_on_interface_spherical(1) = RA;
            this->state_on_interface_spherical(2) = DEC;

            //Step 2: extract mass variable
            this->state_on_interface_spherical(6) = X[Xindex++];
            this->state_on_interface_cartesian(6) = this->state_on_interface_spherical(6);

            //Step 3: left epoch
            this->process_left_epoch(X, Xindex, F, Findex, G, needG);
            this->state_on_interface_spherical(7) = this->EventLeftEpoch;
            this->state_on_interface_cartesian(7) = this->EventLeftEpoch;

            //Step 4: convert to cartesian

            this->state_on_interface_cartesian(0) = r * cosRA * cosDEC;
            this->state_on_interface_cartesian(1) = r * sinRA * cosDEC;
            this->state_on_interface_cartesian(2) = r * sinDEC;

            //Step 5: derivatives
            if (needG)
            {
                //position derivatives
                double r3 = (r * r * r) _GETVALUE;
                double dr_dRA = ( (cosDEC2 * cosRA * sinRA) * (1 / a2 - 1 / b2) )_GETVALUE / r3;
                double dr_dDEC = ( (cosDEC * sinDEC * cosRA2) / a2 + (cosDEC * sinDEC * sinRA2) / b2 - (cosDEC * sinDEC) / c2 )_GETVALUE / r3;
                double dx_dRA = (-r * sinRA * cosDEC + cosRA * cosDEC * dr_dRA)_GETVALUE;
                double dx_dDEC = (-r * cosRA * sinDEC + cosRA * cosDEC * dr_dDEC)_GETVALUE;
                double dy_dRA = (r * cosRA * cosDEC + sinRA * cosDEC * dr_dRA)_GETVALUE;
                double dy_dDEC = (-r * sinRA * sinDEC + sinRA * cosDEC * dr_dDEC)_GETVALUE;
                double dz_dRA = (sinDEC * dr_dRA) _GETVALUE;
                double dz_dDEC = (r * cosDEC + sinDEC * dr_dDEC)_GETVALUE;

                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_RA[0]]) = dx_dRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_RA[1]]) = dy_dRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_RA[2]]) = dz_dRA;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_DEC[0]]) = dx_dDEC;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_DEC[1]]) = dy_dDEC;
                std::get<2>(this->Derivatives_of_state_on_interface_cartesian[this->dIndex_interface_position_wrt_DEC[2]]) = dz_dDEC;
            }
        }//end process_event_right_side

    }//close namespace events
}//close namespace EMTG