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

//EMTGv9 Journey class
//Jacob Englander 6-22-2017


#pragma once

#include "doubleType.h"

#include <vector>

#include "missionoptions.h"
#include "phase.h"
#include "mission.h"
#include "universe.h"
#include "EMTG_enums.h"
#include "LaunchVehicle.h"
#include "Spacecraft.h"

#include "boost/ptr_container/ptr_vector.hpp"

#include "sparsey_thing.h"

namespace EMTG 
{
    //forward declaration of phase, allows "pointer to child" to work
    class phase;

    //forward declaration of mission, allows "pointer to parent" to work
    class Mission;

    class Journey : public sparsey_thing
    {
    public:
        //constructor
        Journey();
        Journey(const size_t& journeyIndex,
            size_t& stageIndex,
            missionoptions* options,
            Astrodynamics::universe* Universe,
            Journey* PreviousJourney,
            Mission* myMission,
            HardwareModels::LaunchVehicle* launchvehicle,
            HardwareModels::Spacecraft* spacecraft);

        //destructor
        virtual ~Journey() {};

        //methods
        //evaluate function
        //return 0 if successful, 1 if failure
        void process_journey( const std::vector<doubleType>& X,
                       size_t& Xindex,
                       std::vector<doubleType>& F,
                       size_t& Findex,
                       std::vector<double>& G,
                       const bool& needG);

        //output functions
        //main output function
        void output(std::ofstream& outputfile,
                    size_t& jprint,
                    size_t& eventcount);

        void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

        void output_STMs();

        void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

        //method to output journey header
        void output_journey_header( std::ofstream& outputfile,
                                    size_t& jprint);

        //bounds calculation function
        //return 0 if successful, 1 if failure
        void calcbounds(std::vector<double>& Xupperbounds,
                        std::vector<double>& Xlowerbounds,
                        std::vector<double>& X_scale_factors,
                        std::vector<double>& Fupperbounds, 
                        std::vector<double>& Flowerbounds,
						std::vector<double>& F_scale_factors,
                        std::vector<std::string>& Xdescriptions,
                        std::vector<std::string>& Fdescriptions,
                        std::vector<size_t>& iGfun,
                        std::vector<size_t>& jGvar,
                        std::vector<std::string>& Gdescriptions,
                        std::vector<size_t>& iAfun,
                        std::vector<size_t>& jAvar,
                        std::vector<std::string>& Adescriptions,
                        std::vector<double>& A,
                        std::vector<double>& synodic_periods);

        //get
        size_t getNumberOfPhases() const { return this->number_of_phases; }
        Phases::phase* getPhase(const size_t& phaseIndex) { return &this->phases[phaseIndex]; }
        Phases::phase* getFirstPhase() { return &this->phases.front(); }
        Phases::phase* getLastPhase() { return &this->phases.back(); }
        doubleType getDeterministicDeltav() const;
        doubleType getStatisticalDeltav() const;
        doubleType getFinalMass() const { return this->phases.back().getFinalMass(); }

    private:
        //pointer to previous journey
        Journey* PreviousJourney;

        //vector of phases
        boost::ptr_vector< Phases::phase > phases;

        //vector of boundary states
        std::vector< math::Matrix<doubleType> > boundary_states;

        //fields
        std::vector<std::string> phase_type_codes;
        size_t number_of_phases;
        size_t journeyIndex;

        //pointer to options
        missionoptions* myOptions;
        JourneyOptions* myJourneyOptions;

        //pointer to parent
        Mission* myMission;

        //pointer to Universe
        Astrodynamics::universe* myUniverse;

        //pointer to hardware
        HardwareModels::LaunchVehicle* myLaunchVehicle;
        HardwareModels::Spacecraft* mySpacecraft;

        //derivative information
        size_t first_X_entry_in_journey;
        size_t first_F_entry_in_journey;
        size_t first_G_entry_in_journey;
        std::vector<size_t> timeconstraints_G_indices;
        std::vector<size_t> departure_date_constraint_G_indices;
    };

} /* namespace EMTG */


