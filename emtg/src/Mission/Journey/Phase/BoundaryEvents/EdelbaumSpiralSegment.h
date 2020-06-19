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

//header for Edelbaum spiral segment
//Jacob Englander 10-24-2017

#pragma once

#include "doubleType.h"

#include "Spacecraft.h"
#include "EMTG_Matrix.h"
#include "universe.h"
#include "missionoptions.h"
#include "EdelbaumSpiral.h"
#include "writey_thing.h"
#include "sparsey_thing.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //forward declare EdelbaumSpiral
        class EdelbaumSpiral;

        class EdelbaumSpiralSegment : public writey_thing, public sparsey_thing
        {
        public:
            EdelbaumSpiralSegment();

            EdelbaumSpiralSegment(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                const size_t& spiralSegmentIndex,
                Astrodynamics::universe* myUniverse,
                Astrodynamics::body* myBody,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                EdelbaumSpiral* mySpiral);

            //output
            void output(std::ofstream& outputfile,
                size_t& eventcount);

            //calcbounds goes in the specialized phase
            void calcbounds();

            void setup_calcbounds(std::vector<double>* Xupperbounds,
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

            //process
            void process(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //get/set
            void set_previous_segment(EdelbaumSpiralSegment* previousSegment) { this->previousSegment = previousSegment; }
            void set_radius_before_segment(double& radius_before_segment) { this->radius_before_segment = &radius_before_segment; }
            void set_radius_after_segment(double& radius_after_segment) { this->radius_after_segment = &radius_after_segment; }
            void set_deltav(const double& deltav) { this->deltav = deltav; }
            doubleType get_mass_after_segment() const { return this->mass_after_segment; }
            size_t get_Xindex_mass_after_segment() const { return this->Xindex_mass_after_segment; }
            size_t get_Xindex_spiralSegmentTime() const { return this->Xindex_SpiralSegmentTime; }

        protected:
            void calculate_dependencies_epoch_time();

            virtual void process_epoch_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            size_t spiralSegmentIndex;
            Astrodynamics::universe* myUniverse;
            Astrodynamics::body* myBody;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;
            JourneyOptions* myJourneyOptions;
            EdelbaumSpiral* mySpiral;
            
            doubleType deltav;
            doubleType LaunchDate;
            doubleType segmentLeftEpoch;
            bool isHeliocentric;

            math::Matrix<doubleType> StateRelativeToJourneyCentralBody;
            math::Matrix<doubleType> PositionRelativeToSun;
            math::Matrix<double> dPositionRelativeToSun_dt;
            math::Matrix<doubleType> PositionCB_Sun;
            math::Matrix<double> dPositionCB_Sun_dt;

            double myDutyCycle;

            doubleType mass_before_segment;
            doubleType mass_after_segment;

            double* radius_before_segment;
            double* radius_after_segment;

            doubleType SpiralSegmentTime;
            doubleType virtual_electric_propellant;
            doubleType virtual_chemical_fuel;

            doubleType actual_SpiralSegmentTime;
            doubleType electric_propellant_used;
            doubleType chemical_fuel_used;

            size_t Xindex_mass_after_segment;
            size_t Xindex_SpiralSegmentTime;
            size_t Xindex_virtual_electric_propellant;
            size_t Xindex_virtual_chemical_fuel;

            size_t Findex_mass_after_segment;
            size_t Findex_SpiralSegmentTime;
            size_t Findex_virtual_electric_propellant;
            size_t Findex_virtual_chemical_fuel;

            //the time continuity constraint has derivatives with respect to all previous time variables, left-hand mass, and segment encoded time
            std::vector<size_t> Gindex_SpiralSegmentTime_wrt_previous_time_variables;
            size_t Gindex_SpiralSegmentTime_wrt_left_hand_mass;
            size_t Gindex_SpiralSegmentTime_wrt_SpiralSegmentTime;

            //the mass continuity constraint has derivatives with respect to all previous time variables, left-hand mass, segment encoded mass, and IF ACS IS ON also segment encoded time
            std::vector<size_t> Gindex_mass_after_segment_wrt_previous_time_variables;
            size_t Gindex_mass_after_segment_wrt_left_hand_mass;
            size_t Gindex_mass_after_segment_wrt_mass_after_segment;
            size_t Gindex_mass_after_segment_wrt_SpiralSegmentTime;

            //the electric propellant continuity constraint has derivatives with respect to all previous time variables, left-hand mass, and segment encoded electric propellant
            std::vector<size_t> Gindex_virtual_electric_propellant_wrt_previous_time_variables;
            size_t Gindex_virtual_electric_propellant_wrt_left_hand_mass;
            size_t Gindex_virtual_electric_propellant_wrt_virtual_electric_propellant;

            //the chemical fuel continuity constraint has derivatives with respect to segment encoded time and segment encoded chemical fuel
            size_t Gindex_virtual_chemical_fuel_wrt_virtual_chemical_fuel;
            size_t Gindex_virtual_chemical_fuel_wrt_SpiralSegmentTime;

            EdelbaumSpiralSegment* previousSegment;
            std::vector< std::size_t > Xindices_EventLeftEpoch;
            size_t Xindex_left_mass;
        };//end class EdelbaumSpiralSegment
    }//close namespace BoundaryEvents
}//close namespace EMTG