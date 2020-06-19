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

//base subphase for EMTGv9 MGAnDSMs
//Jacob Englander 8-23-2017

#pragma once

#include "doubleType.h"

#include "Spacecraft.h"
#include "EMTG_Matrix.h"
#include "universe.h"
#include "missionoptions.h"
#include "PropagatorBase.h"

#include "MGAnDSMs_maneuver_constraint.h"
#include "IntegrationScheme.h"
#include "TimeDomainSpacecraftEOM.h"

#include "sparsey_thing.h"
#include "writey_thing.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare base maneuver epoch constraint object
        class MGAnDSMs_maneuver_constraint;

        //forward declare phase
        class MGAnDSMs_phase;

        class MGAnDSMs_subphase : virtual public writey_thing, virtual public sparsey_thing
        {
        public:
            //constructor
            MGAnDSMs_subphase();
            MGAnDSMs_subphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& subphaseIndex,
                const size_t& stageIndex,
                MGAnDSMs_phase* myPhase,
                MGAnDSMs_subphase* previousSubPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual MGAnDSMs_subphase* clone() const = 0;

            //propagator
            virtual ~MGAnDSMs_subphase();

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount) = 0;

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream & acceleration_model_file) = 0;

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget) = 0;

            //calcbounds goes in the specialized phase
            virtual void calcbounds(std::vector<size_t>& timeVariables) = 0;

            void calcbounds_maneuver_constraints();

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

            //process goes in the specialized phase
            virtual void process_subphase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            void process_maneuver_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //configure propagator
            virtual void configure_propagator() = 0;

            //get/set
            inline void set_spacecraft_state_minus(math::Matrix<doubleType>& spacecraft_state_minus) { this->spacecraft_state_minus_pointer = &spacecraft_state_minus; }
            inline void set_spacecraft_state_plus(math::Matrix<doubleType>& spacecraft_state_plus) { this->spacecraft_state_plus_pointer = &spacecraft_state_plus; }
            inline void set_SPTM(math::Matrix<double>& SPTM) { this->SPTM_pointer = &SPTM; }
            inline void set_dPropagatedStatedIndependentVariable(math::Matrix<double>& dPropagatedStatedIndependentVariable_pointer) { this->dPropagatedStatedIndependentVariable_pointer = &dPropagatedStatedIndependentVariable_pointer; }
            inline void set_PhaseFlightTime(doubleType& PhaseFlightTime) { this->PhaseFlightTime_pointer = &PhaseFlightTime; }
            inline void set_timeVariables(std::vector<size_t>& timeVariables) { this->timeVariables = timeVariables; }
            inline doubleType getDSMmagnitude() const { return this->DSM_magnitude; }
            inline doubleType getTCMmagnitude() const { return this->TCM_magnitude; }
            inline std::string getName() const { return this->name; }
            inline std::vector<size_t>& getXindex_DSM_components() { return this->Xindex_DSM_components; }
            inline math::Matrix<doubleType>& getDSM() { return this->DSM; }
            inline std::vector<size_t>& getXindex_burnIndex() { return this->Xindex_BurnIndex; }
            inline doubleType getDSMepoch() const { return this->DSMepoch; }
            inline doubleType getSubPhaseTime() const { return this->SubPhaseTime; }
            inline size_t getFirst_X_entry_this_subphase() const { return this->First_X_entry_this_subphase; }
            inline bool getBordersMatchPoint() const { return this->BordersMatchPoint; }
            virtual bool getBordersBoundary() const = 0;
            inline std::vector<size_t> get_timeVariables() { return this->timeVariables; }
            inline MGAnDSMs_subphase* get_previousSubphase() { return this->previousSubPhase; }
            inline MGAnDSMs_phase* get_myPhase() { return this->myPhase; }

            inline math::Matrix<double> get_dFuel_dDSMcomponents() const { return this->dFuel_dDSMcomponents; }
            inline math::Matrix<double> get_dOxidizer_dDSMcomponents() const { return this->dOxidizer_dDSMcomponents; }

        protected:
            //calcbounds
            virtual void calcbounds_subphase_time() = 0;

            void calcbounds_DSM_components();

            //process

            virtual void process_subphase_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            void process_DSM_components(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void write_output_line(std::ofstream& outputfile,
                size_t& eventcount,
                const std::string& event_type,
                const std::string& location,
                const doubleType& timestep_size,
                const doubleType& angle1,
                const doubleType& angle2,
                const math::Matrix<doubleType>& state,
                const math::Matrix<doubleType>& dV,
                const doubleType dVmag,
                const doubleType Isp); //override what is in writey_thing() because we don't need power/propulsion stuff

            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            size_t subphaseIndex;
            MGAnDSMs_subphase* previousSubPhase;
            Astrodynamics::universe* myUniverse;
            HardwareModels::Spacecraft* mySpacecraft;
            JourneyOptions* myJourneyOptions;
            MGAnDSMs_phase* myPhase;
            bool BordersMatchPoint;
            PropagatorType myPropagatorType;

            std::string prefix;

            size_t First_X_entry_this_subphase;

            //pointers to external things
            math::Matrix<doubleType>* spacecraft_state_minus_pointer;
            math::Matrix<doubleType>* spacecraft_state_plus_pointer;
            doubleType* PhaseFlightTime_pointer;
            math::Matrix<double>* SPTM_pointer;
            
            //internal states
            math::Matrix<doubleType> StateAfterPropagationBeforeDSM;
            math::Matrix<doubleType> StateAfterDSMBeforeTCM;

            //time
            std::vector<size_t> timeVariables;

            //propagator
            bool isKeplerian;
            Astrodynamics::PropagatorBase* myPropagator;
            math::Matrix<double> STM;
            math::Matrix<double>* dPropagatedStatedIndependentVariable_pointer;
            doubleType BurnIndex;
            double BurnIndexDouble;
            std::vector<size_t> Xindex_BurnIndex;//last entry is the one we care about
            doubleType SubPhaseTime;            
            
            //maneuvers
            math::Matrix<double> MTM;
            std::vector<size_t> Xindex_DSM_components;
            math::Matrix<doubleType> DSM;
            std::vector<size_t> Gindices_DSM_magnitude_constraint;
            doubleType TCM_magnitude;
            doubleType DSM_magnitude;
            doubleType DSMepoch;
            bool hasTCM;
            boost::ptr_vector<MGAnDSMs_maneuver_constraint> myManeuverConstraints;
            PropulsionSystemChoice ChemicalManeuverType;

            //tanks
            doubleType chemical_fuel_used;
            doubleType chemical_oxidizer_used;
            doubleType DSM_fuel_used;
            doubleType TCM_fuel_used;
            doubleType ACS_fuel_used;
            math::Matrix<double> dFuel_dDSMcomponents;
            math::Matrix<double> dOxidizer_dDSMcomponents;

            //integrator pointers
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;

            //ephemeris resolution
            double EphemerisOutputResolution;

            //output things
            size_t num_timesteps;
            math::Matrix<doubleType> output_state;
            math::Matrix<doubleType> temp_state;
            math::Matrix<doubleType> empty3;
        };

        inline MGAnDSMs_subphase * new_clone(MGAnDSMs_subphase const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG