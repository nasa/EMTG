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

#include "EphemerisPeggedFlybyOut.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedPoweredFlyby : public EphemerisPeggedFlybyOut
        {
        public:
            //specialized constructor
            EphemerisPeggedPoweredFlyby(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
				EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent);

            //output
            void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

            math::Matrix<doubleType> get_periapse_state();

            //process
            void process_event(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //calcbounds
            void calcbounds(std::vector<size_t> timeVariables);

            //output
            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

        private:
            void calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds);

            void calcbounds_virtual_propellant_constraints();

            void calcbounds_deltav_contribution();

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            size_t dIndex_Mass_wrt_FlybyPeriapseDistance;
            std::vector<size_t> dIndices_Mass_wrt_Vinfinity_out;
            std::vector<size_t> dIndices_Mass_wrt_Vinfinity_in;

            doubleType FlybyDeltavSigned;
            doubleType FlybyDeltavMagnitude;
            doubleType FlybyPeriapseDistance;
            
            double dChemicalFuel_dDeltavMagnitude;
            double dChemicalOxidizer_dDeltavMagnitude;
            math::Matrix<double> dDeltavMagnitude_dVinfinity_in;
            math::Matrix<double> dDeltavMagnitude_dVinfinity_out;
            double dDeltavMagnitude_dFlybyPeriapseDistance;
            double dChemicalFuel_dInitialMass;
            double dChemicalOxidizer_dInitialMass;

            std::vector< std::vector<size_t> > Gindices_TurnAngleConstraint_with_respect_to_Vinfinity; //first incoming, then outgoing
            size_t Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance;

            std::vector< std::vector<size_t> > Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity; //first incoming, then outgoing
            size_t Gindices_VirtualChemicalFuel_with_respect_to_FlybyPeriapseDistance;
            size_t Gindices_dVirtualChemicalFuel_dLeftMass;

            std::vector< std::vector<size_t> > Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity; //first incoming, then outgoing
            size_t Gindices_VirtualChemicalOxidizer_with_respect_to_FlybyPeriapseDistance;
            size_t Gindices_dVirtualChemicalOxidizer_dLeftMass;

            std::vector< std::vector<size_t> > Gindices_Deltav_with_respect_to_Vinfinity; //first incoming, then outgoing
            size_t Gindices_Deltav_with_respect_to_FlybyPeriapseDistance;
        };
    }//end namespace BoundaryEvents
}//end namespace EMTG