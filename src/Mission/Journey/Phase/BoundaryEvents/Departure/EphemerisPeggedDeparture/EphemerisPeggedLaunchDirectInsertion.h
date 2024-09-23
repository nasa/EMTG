
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

#include "EphemerisPeggedDeparture.h"
#include "LaunchVehicle.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedLaunchDirectInsertion : public EphemerisPeggedDeparture
        {
        public:
            //specialized constructor
            EphemerisPeggedLaunchDirectInsertion(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions,
                ArrivalEvent* PreviousPhaseArrivalEvent);

            //output
            void output(std::ofstream& outputfile,
                        const double& launchdate,
                        size_t& eventcount);

            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            math::Matrix<doubleType> get_periapse_state() { return math::Matrix<doubleType>(6, 1, 0.0); };//TODO: eventually I might use this to derive a periapse state?

            inline size_t getXindex_of_initial_impulse() const { return this->Xindex_C3; }

            //process
            void process_event(const std::vector<doubleType>& X,
                               size_t& Xindex,
                               std::vector<doubleType>& F,
                               size_t& Findex,
                               std::vector<double>& G,
                               const bool& needG);

            //calcbounds
            void calcbounds(std::vector<size_t> timeVariables);

        private:
            void calcbounds_event_left_side(std::vector<size_t> timeVariables);

            void calcbounds_event_main() {}; //stub - all of the interesting stuff happens in calcbounds_event_left_side()

            void calcbounds_virtual_propellant_constraints(); 
            
            void calcbounds_deltav_contribution();

            void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG); //stub - all of the interesting stuff happens in process_event_left_side()

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            EMTG::HardwareModels::LaunchVehicle* myLaunchVehicle;


            //fields
            bool useLV;
            bool encodeInitialMass;
            doubleType RLA;
            doubleType DLA;
            size_t Xindex_InitialMassMultiplier;
            size_t Xindex_C3;
            size_t Xindex_RLA;
            size_t Xindex_DLA;
            size_t dIndex_C3_vx;
            size_t dIndex_C3_vy;
            size_t dIndex_C3_vz;
            size_t dIndex_C3_mass;
            size_t dIndex_RLA_vx;
            size_t dIndex_RLA_vy;
            size_t dIndex_DLA_vx;
            size_t dIndex_DLA_vy;
            size_t dIndex_DLA_vz;            
            size_t dIndex_InitialMassMultiplier_mass;

            //derivative entries
            double dm_dVinfinityOut;
            doubleType InitialMassMultiplier;
            doubleType MassBeforeMultiplier;

            //delta-v
            size_t Gindex_dDeltav_dVinfinity;
            
            //propellant
            size_t Gindices_dVirtualChemicalFuel_dLeftMass;
            size_t Gindices_dVirtualChemicalFuel_dVinfinity;
            size_t Gindices_dVirtualChemicalOxidizer_dLeftMass;
            size_t Gindices_dVirtualChemicalOxidizer_dVinfinity;

            double dChemicalFuel_dVinfinity;
            double dChemicalOxidizer_dVinfinity;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG